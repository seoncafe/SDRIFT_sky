#!/usr/bin/env python3
"""Query and plot all S-DRIFT sky models through one API using the docs/data files."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import warnings
from pathlib import Path
import numpy as np
import healpy as hp
from sky_coordinates import COORD_SYSTEMS, check_coords, from_equatorial, to_equatorial

DEFAULT_DATA = Path(__file__).resolve().parent / 'docs' / 'data'
BACKGROUND_MODELS = ('zodiacal', 'dgl', 'zodiacal_dgl')
MODELS = BACKGROUND_MODELS + ('airglow', 'total')
BANDS = ('UV1', 'UV2', 'UV3', 'VIS')
ARCSEC_PER_RADIAN = 206264.80624709636


def utc_time(date):
    stamp = datetime.fromisoformat(str(date).replace('Z', '+00:00'))
    if stamp.tzinfo is None:
        stamp = stamp.replace(tzinfo=timezone.utc)
    stamp = stamp.astimezone(timezone.utc)
    return stamp.timestamp(), stamp.isoformat(timespec='seconds').replace('+00:00', 'Z')


def vector(ra, dec):
    if not np.isfinite(ra) or not np.isfinite(dec) or not -90 <= dec <= 90:
        raise ValueError('RA must be finite; DEC must be between -90 and 90 degrees')
    a, d = np.deg2rad([ra % 360, dec])
    return np.array([np.cos(d)*np.cos(a), np.cos(d)*np.sin(a), np.sin(d)])


def coordinates(v):
    return float(np.rad2deg(np.arctan2(v[1], v[0])) % 360), float(np.rad2deg(np.arcsin(np.clip(v[2], -1, 1))))


def rotation(v, angle):
    # Same fixed J2000 ecliptic rotation as docs/assets/sky-core.mjs.
    c, s = np.cos(np.deg2rad(23.439291111)), np.sin(np.deg2rad(23.439291111))
    x, y, z = np.moveaxis(np.asarray(v), -1, 0)
    ye, ze = c*y+s*z, -s*y+c*z
    xr, yr = np.cos(angle)*x-np.sin(angle)*ye, np.sin(angle)*x+np.cos(angle)*ye
    return np.stack([xr, c*yr-s*ze, s*yr+c*ze], axis=-1)


def sun_longitude(v):
    e = np.deg2rad(23.439291111)
    return np.arctan2(np.cos(e)*v[1]+np.sin(e)*v[2], v[0])


def zodiacal_template(coeffs, sun_lon):
    """Combine annual harmonics (A0, C1, S1, ...) at Sun J2000 ecliptic longitude.

    Same as zodiacal_healpix.harmonic_weights and docs/assets/sky-core.mjs.
    """
    coeffs = np.asarray(coeffs, dtype=float)
    weights = [1.]
    for k in range(1, (len(coeffs)-1)//2+1):
        weights += [np.cos(k*sun_lon), np.sin(k*sun_lon)]
    return np.asarray(weights)@coeffs


def ab_magnitude(flux):
    with np.errstate(divide='ignore', invalid='ignore'):
        return -2.5*np.log10(np.asarray(flux)*1e6/ARCSEC_PER_RADIAN**2/3631.)


def boundary(v, angle):
    ref = np.array([0., 0., 1.]) if abs(v[2]) < .9 else np.array([1., 0., 0.])
    a = np.cross(v, ref); a /= np.linalg.norm(a)
    b = np.cross(v, a)
    t = np.linspace(0, 2*np.pi, 721); r = np.deg2rad(angle)
    return np.cos(r)*v[:, None] + np.sin(r)*(a[:, None]*np.cos(t)+b[:, None]*np.sin(t))


class _BackgroundModel:
    """Cached model data; coordinates are equatorial degrees in fixed GCRS axes."""
    def __init__(self, data_dir=DEFAULT_DATA):
        self.data_dir = Path(data_dir)
        self.manifest = json.loads((self.data_dir/'manifest.json').read_text())
        self._cache = {}

    def _array(self, meta):
        name = meta['file']
        if name not in self._cache:
            raw = gzip.decompress((self.data_dir/name).read_bytes())
            if len(raw) != meta['bytes'] or hashlib.sha256(raw).hexdigest() != meta['sha256']:
                raise ValueError(f'Dataset size/SHA-256 mismatch: {name}')
            self._cache[name] = np.frombuffer(raw, dtype='<f4').reshape(meta['shape'])
        return self._cache[name]

    def _body_vector(self, unix, name):
        """Linearly interpolated geocentric vector (unit for the Sun, km for the Moon)."""
        meta = self.manifest[name]
        u = (unix-meta['start_unix'])/meta['step_seconds']; i = int(np.floor(u)); f = u-i
        if i < 0 or i >= meta['count']-1:
            raise ValueError('UTC date must be within 1900-2100')
        data = self._array(meta)
        return data[i].astype(float)*(1-f)+data[i+1].astype(float)*f

    def _body(self, unix, name):
        v = self._body_vector(unix, name)
        return v/np.linalg.norm(v)

    def evaluate(self, directions, date, band='VIS', model='zodiacal_dgl',
                 sun_avoidance=True, sun_separation=90., moon_avoidance=True, moon_separation=60.):
        """Return unmasked component fluxes plus a common avoidance validity mask."""
        if model not in BACKGROUND_MODELS or band not in BANDS:
            raise ValueError('Unknown model or band')
        if not all(np.isfinite(a) and 0 <= a <= 180 for a in (sun_separation, moon_separation)):
            raise ValueError('Separation angles must be finite and within 0–180 degrees')
        v = np.asarray(directions, dtype=float)
        if v.shape[-1] != 3 or not np.all(np.isfinite(v)) or not np.allclose(np.linalg.norm(v, axis=-1), 1):
            raise ValueError('Directions must be finite unit vectors with shape (..., 3)')
        unix, utc = utc_time(date)
        sun, moon = self._body(unix, 'solar'), self._body(unix, 'lunar')
        sun_lon = sun_longitude(sun)
        delta = (sun_lon-sun_longitude(self.manifest['solar']['reference_vector'])+np.pi)%(2*np.pi)-np.pi
        def sample(data, direction):
            return hp.get_interp_val(data, np.arccos(np.clip(direction[..., 2], -1, 1)),
                                     np.mod(np.arctan2(direction[..., 1], direction[..., 0]), 2*np.pi))
        if model != 'dgl':
            template = zodiacal_template(self._array(self.manifest['zodiacal'][band]), sun_lon)
            zodi = sample(template, rotation(v, -delta))
        else:
            zodi = np.zeros(v.shape[:-1])
        dgl = sample(self._array(self.manifest['dgl']['bands'][band]), v) if model != 'zodiacal' else np.zeros(v.shape[:-1])
        sun_dot, moon_dot = v@sun, v@moon
        valid = ((sun_dot <= np.cos(np.deg2rad(sun_separation))) | (not sun_avoidance)) & ((moon_dot <= np.cos(np.deg2rad(moon_separation))) | (not moon_avoidance))
        flux = zodi+dgl
        return dict(utc=utc, band=band, model=model, zodiacal=zodi, dgl=dgl, flux=flux,
                    mag=ab_magnitude(flux), valid=valid, sun=sun, moon=moon, delta=delta if model != 'dgl' else 0.,
                    sun_angle=np.rad2deg(np.arccos(np.clip(sun_dot, -1, 1))),
                    moon_angle=np.rad2deg(np.arccos(np.clip(moon_dot, -1, 1))),
                    sun_avoidance=bool(sun_avoidance), sun_separation=float(sun_separation),
                    moon_avoidance=bool(moon_avoidance), moon_separation=float(moon_separation))

    def point(self, ra, dec, date, **settings):
        """Direct interpolated sightline query, including unmasked brightness."""
        r = self.evaluate(vector(ra, dec), date, **settings)
        visible = bool(r['valid']); mag = float(r['mag'])
        return dict(ra_deg=float(ra % 360), dec_deg=float(dec), utc=r['utc'], band=r['band'], model=r['model'],
                    ab_mag_arcsec2=mag, intensity_mjy_sr=float(r['flux']),
                    zodiacal_mjy_sr=float(r['zodiacal']), dgl_mjy_sr=float(r['dgl']),
                    visible=visible, display_ab_mag_arcsec2=mag if visible else None,
                    sun_ra_dec_deg=coordinates(r['sun']), moon_ra_dec_deg=coordinates(r['moon']),
                    sun_separation_deg=float(r['sun_angle']), moon_separation_deg=float(r['moon_angle']),
                    sun_avoidance=r['sun_avoidance'], minimum_sun_separation_deg=r['sun_separation'],
                    moon_avoidance=r['moon_avoidance'], minimum_moon_separation_deg=r['moon_separation'])

    def sky_map(self, date, nside=128, coords='equatorial', **settings):
        check_coords(coords)
        if not hp.isnsideok(nside, nest=True):
            raise ValueError('NSIDE must be a positive power of two')
        xyz = np.stack(hp.pix2vec(nside, np.arange(hp.nside2npix(nside))), axis=-1)
        r = self.evaluate(to_equatorial(xyz, coords), date, **settings)
        r['nside'] = nside
        r['coords'] = coords
        return r

    def save_fits(self, result, filename):
        r = result
        sun_ra, sun_dec = coordinates(r['sun']); moon_ra, moon_dec = coordinates(r['moon'])
        arrays = [np.where(r['valid'], r['mag'], np.nan), np.where(r['valid'], r['flux'], np.nan)]
        names = ['AB_MAG_ARCSEC2', 'INTENSITY_MJY_SR']; units = ['mag/arcsec2', 'MJy/sr']
        if r['model'] == 'zodiacal_dgl':
            arrays += [np.where(r['valid'], r[k], np.nan) for k in ('zodiacal', 'dgl')]
            names += ['ZODIACAL_MJY_SR', 'DGL_MJY_SR']; units += ['MJy/sr', 'MJy/sr']
        lo, hi = self.manifest['bands'][r['band']]
        headers = [('DATE-UTC', r['utc']), ('BAND', r['band']), ('BANDMIN', lo/1000), ('BANDMAX', hi/1000),
                   ('DATATYPE', r['model']), ('ANTISUN', r['sun_avoidance']), ('ELONGMIN', r['sun_separation']),
                   ('SUNRA', sun_ra), ('SUNDEC', sun_dec), ('MOONMASK', r['moon_avoidance']),
                   ('MOONMIN', r['moon_separation']), ('MOONRA', moon_ra), ('MOONDEC', moon_dec),
                   ('MOONOBS', 'SATELLITE' if 'state' in r else 'GEOCENTRIC'), ('ROT-DEG', np.rad2deg(r['delta'])), ('PHOTCNT', True)]
        if r['model'] != 'dgl': headers += [('REF-UTC', self.manifest['solar']['reference_utc'])]
        if r['model'] != 'zodiacal':
            meta = self.manifest['dgl']['metadata']
            headers += [('CORRMOD', 'C2022_all'), ('CORRSCL', meta['scale']), ('BETAOP', 'MULTIPLY'),
                        ('ALPHAPP', meta['alpha_approximation']), ('BETAREF', meta['beta_reference']), ('DGLFIXED', True)]
        if 'state' in r:
            keys=['earth','valid','hmin','theta','view','elv']
            arrays += [r[k] for k in keys]
            names += ['EARTH_OCCULTED','MODEL_VALID','MIN_ALTITUDE','ZENITH_ANGLE','VIEW_ALLOWED','ELV']
            units += ['','','km','deg','','deg']
            state=r['state']; zen_ra,zen_dec=coordinates(state['zenith'])
            headers += [('ORBIT',r['orbit']),('ORBMODEL','SSO J2 CIRCULAR'),('SAT_H',r['height']),
                        ('INCL',state['inclination_deg']),('RAAN',state['raan_deg']),('LTAN',state['ltan_hours']),
                        ('SAT_SZA',state['sza']),('PHASE',state['orbital_phase_deg']),('PHASE0',r['orbital_phase_deg']),
                        ('PHREF','NORTH'),('PHCONV','NORTH_V2'),('PHMODE','EPOCH' if r['epoch'] else 'MANUAL'),
                        ('PERIOD',state['period']),('ELVMIN',r['min_elv_deg']),('ELVREF','SPHERICAL EARTH'),
                        ('EARTHR',6371.),('ZENRA',zen_ra),('ZENDEC',zen_dec)]
            if r['epoch']:headers.append(('EPOCH',r['epoch']))
        path = Path(filename)
        if not path.name.endswith('.fits.gz'): raise ValueError('FITS output must end with .fits.gz')
        path.parent.mkdir(parents=True, exist_ok=True)
        hp.write_map(str(path), arrays, column_names=names, column_units=units, dtype=np.float32,
                     coord='G' if r.get('coords')=='galactic' else 'C', overwrite=True, extra_header=headers)


class SkyModel:
    """Common API: all five models, orbit masks, and MJy/sr intensities."""
    def __init__(self, data_dir=DEFAULT_DATA):
        # Lazy import keeps the numerical engine independent of this facade.
        from airglow_orbit import _OrbitModel
        self._engine = _OrbitModel(data_dir)
        self.data_dir = self._engine.data_dir
        self.manifest = self._engine.manifest

    def evaluate(self, directions, date, model='total', airglow_model='combined', **settings):
        """Evaluate equatorial unit vectors (N,3); return masked common fields.

        model selects sky components; airglow_model selects combined/glow.
        Legacy kind= and model=combined/glow are accepted as aliases.
        """
        legacy_kind = settings.pop('kind', None)
        if model in ('combined', 'glow'):
            airglow_model, model = model, legacy_kind or 'total'
        elif legacy_kind is not None:
            if model != 'total' and model != legacy_kind:
                raise ValueError('model and legacy kind specify different sky components')
            model = legacy_kind
        coords = check_coords(settings.pop('coords', 'equatorial'))
        v = np.asarray(directions, dtype=float)
        if v.ndim not in (1, 2) or v.shape[-1] != 3:
            raise ValueError('Directions must have shape (3,) or (N, 3)')
        v = v.reshape(-1, 3)
        raw = self._engine.evaluate(v, date, kind=model, model=airglow_model, **settings)
        raw['coords'] = coords
        state = raw['state']; height = raw['height']
        b = (6371.+height)*(v@state['zenith'])
        earth = (b < 0) & (b*b-(6371.+height)**2+6371.**2 >= 0)
        hmin = np.sqrt(np.maximum((6371.+height)**2-np.minimum(b, 0)**2, 0))-6371.
        visible = raw['view'] & ~earth
        n = len(v); zeros = np.zeros(n)
        photon = raw['flux'] if model == 'airglow' else raw.get('airglow_photon_radiance', zeros)
        lo, hi = self.manifest['bands'][raw['band']]
        air = 6.62607015e-10*photon/np.log(hi/lo)
        components = {key: np.where(visible, value, np.nan) for key, value in
                      [('zodiacal', raw.get('zodiacal', zeros)), ('dgl', raw.get('dgl', zeros)), ('airglow', air)]}
        flux = components['zodiacal']+components['dgl']+components['airglow']
        result = dict(raw)
        result.update(components)
        result.update(model=model, kind=model, airglow_model=airglow_model, coords=coords,
            flux=flux, intensity_mjy_sr=flux, flux_unit='MJy/sr', mag=ab_magnitude(flux),
            airglow_photon_radiance=np.where(visible, photon, np.nan),
            lit=np.where(visible, raw.get('lit', zeros), np.nan),
            earth=earth, earth_occulted=earth, hmin=np.where(earth, 0., hmin),
            visible=visible, model_valid=np.asarray(raw['valid'], dtype=bool),
            valid=visible & raw['valid'], view_allowed=raw['view'],
            cutoff=np.asarray(raw.get('cutoff', np.zeros(n, bool)), dtype=bool),
            _engine_result=raw)
        result['ab_mag_arcsec2'] = result['mag']
        return result

    def point(self, ra, dec, date, model='total', airglow_model='combined', **settings):
        """A single RA/DEC direction in degrees; the input is always equatorial."""
        r = self.evaluate(vector(ra, dec), date, model=model, airglow_model=airglow_model, **settings)
        sun_angle = float(np.rad2deg(np.arccos(np.clip(vector(ra, dec)@r['sun'], -1, 1))))
        moon_angle = float(np.rad2deg(np.arccos(np.clip(vector(ra, dec)@r['moon'], -1, 1))))
        return dict(ra_deg=float(ra % 360), dec_deg=float(dec), utc=r['utc'], band=r['band'],
            model=r['model'], airglow_model=r['airglow_model'], coords=r['coords'],
            ab_mag_arcsec2=float(r['mag'][0]), intensity_mjy_sr=float(r['flux'][0]),
            zodiacal_mjy_sr=float(r['zodiacal'][0]), dgl_mjy_sr=float(r['dgl'][0]),
            airglow_mjy_sr=float(r['airglow'][0]), photon_radiance=float(r['airglow_photon_radiance'][0]),
            visible=bool(r['visible'][0]), view_allowed=bool(r['view'][0]),
            earth_occulted=bool(r['earth'][0]), model_valid=bool(r['model_valid'][0]), valid=bool(r['valid'][0]),
            display_ab_mag_arcsec2=float(r['mag'][0]) if r['visible'][0] else None,
            orbit=r['orbit'], height_km=float(r['height']), orbital_phase_deg=r['state']['orbital_phase_deg'],
            sza_deg=r['state']['sza'], zenith_deg=float(r['theta'][0]), elv_deg=float(r['elv'][0]),
            minimum_elv_deg=float(r['min_elv_deg']), sun_ra_dec_deg=coordinates(r['sun']),
            moon_ra_dec_deg=coordinates(r['moon']), sun_separation_deg=sun_angle, moon_separation_deg=moon_angle,
            sun_avoidance=bool(r['sun_avoidance']), minimum_sun_separation_deg=float(r['sun_separation']),
            moon_avoidance=bool(r['moon_avoidance']), minimum_moon_separation_deg=float(r['moon_separation']))

    def sky_map(self, date, nside=128, coords='equatorial', model='total', airglow_model='combined', **settings):
        check_coords(coords)
        if not hp.isnsideok(nside, nest=True): raise ValueError('NSIDE must be a positive power of two')
        v = np.stack(hp.pix2vec(nside, np.arange(hp.nside2npix(nside))), axis=-1)
        r = self.evaluate(to_equatorial(v, coords), date, model=model, airglow_model=airglow_model,
                          coords=coords, **settings)
        r['nside'] = nside
        return r

    def plot(self, result, target=None, contour_step=.2, figure=None):
        """Return a Matplotlib Figure in the result's coordinate system."""
        fig = self._engine.plot(result['_engine_result'] | {'nside': result['nside']},
                                contour_step=contour_step, figure=figure)
        if target is not None:
            from_equ = from_equatorial(vector(*target), result['coords'])
            lon, lat = coordinates(from_equ)
            fig.axes[0].plot(-np.deg2rad((lon+180)%360-180), np.deg2rad(lat),
                             marker='x', color='#009cff', markersize=9, markeredgewidth=2)
        return fig

    def save_fits(self, result, filename):
        """Same intensity units and column layout for every model."""
        from astropy.io import fits
        import tempfile
        path = Path(filename)
        if not path.name.endswith('.fits.gz'): raise ValueError('FITS output must end with .fits.gz')
        raw = result['_engine_result'] | {'nside': result['nside']}
        # Preserve source hashes and orbit metadata from the numerical engine.
        with tempfile.TemporaryDirectory() as tmp:
            legacy = Path(tmp)/'metadata.fits.gz'
            self._engine.save_fits(raw, legacy)
            with fits.open(legacy) as hd:
                metadata = hd[1].header.copy()
        keys = ['mag','flux','zodiacal','dgl','airglow','airglow_photon_radiance',
                'earth','model_valid','hmin','theta','view','elv','visible','valid']
        names = ['AB_MAG_ARCSEC2','INTENSITY_MJY_SR','ZODIACAL_MJY_SR','DGL_MJY_SR','AIRGLOW_MJY_SR',
                 'AIRGLOW_PHOTON_RADIANCE','EARTH_OCCULTED','MODEL_VALID','MIN_ALTITUDE',
                 'ZENITH_ANGLE','VIEW_ALLOWED','ELV','VISIBLE','VALID']
        units = ['mag/arcsec2','MJy/sr','MJy/sr','MJy/sr','MJy/sr','photons/s/cm2/sr',
                 '1','1','km','deg','1','deg','1','1']
        columns = [fits.Column(name=name, format='E', unit=unit,
                   array=np.asarray(result[key], dtype=np.float32)) for key,name,unit in zip(keys,names,units)]
        table = fits.BinTableHDU.from_columns(columns, header=metadata)
        table.header['MODEL'] = result['model']; table.header['DATATYPE'] = result['model'].upper()
        table.header['AGMODEL'] = result['airglow_model']; table.header['APIVER'] = 2
        table.header['FLUXUNIT'] = 'MJy/sr'
        path.parent.mkdir(parents=True, exist_ok=True)
        fits.HDUList([fits.PrimaryHDU(),table]).writeto(path, overwrite=True)


def main():
    """One CLI shared by sdrift_sky.py and airglow_orbit.py."""
    parser = argparse.ArgumentParser(description='S-DRIFT: common sky brightness interface')
    parser.add_argument('--gui', action='store_true')
    parser.add_argument('--date', '--utc', default=datetime.now(timezone.utc).isoformat())
    parser.add_argument('--model', default=None, help='zodiacal, dgl, zodiacal_dgl, airglow, total (default: total)')
    parser.add_argument('--kind', help='Compatibility alias for --model')
    parser.add_argument('--airglow-model', choices=['combined','glow'], default='combined')
    parser.add_argument('--orbit', choices=['dawn-dusk','noon-midnight'], default='dawn-dusk')
    parser.add_argument('--height', '--H', type=float, default=650.)
    parser.add_argument('--orbital-phase-deg', '--phase', type=float, default=0.)
    parser.add_argument('--epoch')
    parser.add_argument('--min-elv', '--elv-min', type=float, default=20.)
    parser.add_argument('--band', choices=BANDS, default='VIS')
    parser.add_argument('--coords', choices=COORD_SYSTEMS, default='equatorial')
    parser.add_argument('--nside', type=int, default=128)
    parser.add_argument('--data-dir', type=Path, default=DEFAULT_DATA)
    parser.add_argument('--engine', choices=['auto','numba','numpy'], default='auto')
    parser.add_argument('--sun-separation', type=float, default=90.)
    parser.add_argument('--moon-separation', type=float, default=60.)
    parser.add_argument('--no-sun-avoidance', action='store_true')
    parser.add_argument('--no-moon-avoidance', action='store_true')
    parser.add_argument('--contour-step', type=float, choices=[i/10 for i in range(1,11)], default=.2)
    parser.add_argument('--ra', type=float); parser.add_argument('--dec', type=float)
    parser.add_argument('--show', action='store_true')
    for name in ['fits','pdf','png','json']: parser.add_argument('--'+name, type=Path)
    args = parser.parse_args()
    try:
        if args.model in ('glow','combined'):
            model, air_model = args.kind or 'total', args.model
        else:
            if args.model and args.kind and args.model != args.kind: raise ValueError('--model and --kind conflict')
            model, air_model = args.model or args.kind or 'total', args.airglow_model
        if model not in MODELS: raise ValueError('Unknown sky model')
        options = dict(model=model, airglow_model=air_model, band=args.band, orbit=args.orbit, height=args.height,
            orbital_phase_deg=args.orbital_phase_deg, epoch=args.epoch, min_elv_deg=args.min_elv,
            sun_avoidance=not args.no_sun_avoidance, sun_separation=args.sun_separation,
            moon_avoidance=not args.no_moon_avoidance, moon_separation=args.moon_separation, engine=args.engine)
        if args.gui:
            from airglow_orbit_gui import main as gui
            gui(args.data_dir, model=model, coords=args.coords, initial=options | {'date':args.date, 'nside':args.nside,
                 'contour_step':args.contour_step})
            return
        sky = SkyModel(args.data_dir)
        if (args.ra is None) != (args.dec is None): raise ValueError('Provide both --ra and --dec')
        if args.ra is not None:
            point = sky.point(args.ra,args.dec,args.date,coords=args.coords,**options)
            print(f"{point['ab_mag_arcsec2']:.6f} AB mag / arcsec² at (RA, DEC) = ({args.ra:g}, {args.dec:g})")
            print(f"model={model}; intensity={point['intensity_mjy_sr']:.8g} MJy/sr; visible={point['visible']}; model_valid={point['model_valid']}")
            if args.json:
                clean={k:None if isinstance(v,float) and not np.isfinite(v) else v for k,v in point.items()}
                args.json.parent.mkdir(parents=True,exist_ok=True)
                args.json.write_text(json.dumps(clean,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
        elif args.json: raise ValueError('--json requires --ra and --dec')
        if args.ra is None or args.show or args.fits or args.png or args.pdf:
            r=sky.sky_map(args.date,nside=args.nside,coords=args.coords,**options)
            print(f"{model} | {r['coords']} | {r['nside']} | visible={np.count_nonzero(r['visible'])} | {r['engine']}")
            if args.fits: sky.save_fits(r,args.fits)
            if args.pdf or args.png or args.show:
                if not args.show:
                    import matplotlib
                    matplotlib.use('Agg')
                fig=sky.plot(r,target=(args.ra,args.dec) if args.ra is not None else None,contour_step=args.contour_step)
                for path in [args.pdf,args.png]:
                    if path:path.parent.mkdir(parents=True,exist_ok=True);fig.savefig(path,dpi=160)
                import matplotlib.pyplot as plt
                if args.show:plt.show()
                else:plt.close(fig)
    except (ValueError,OSError,TypeError,KeyError) as exc: parser.error(str(exc))


if __name__ == '__main__': main()

#!/usr/bin/env python3
"""Query and plot S-DRIFT zodiacal/DGL sky brightness using bundled web data."""
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
MODELS = ('zodiacal', 'dgl', 'zodiacal_dgl')
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


def ab_magnitude(flux):
    with np.errstate(divide='ignore', invalid='ignore'):
        return -2.5*np.log10(np.asarray(flux)*1e6/ARCSEC_PER_RADIAN**2/3631.)


def boundary(v, angle):
    ref = np.array([0., 0., 1.]) if abs(v[2]) < .9 else np.array([1., 0., 0.])
    a = np.cross(v, ref); a /= np.linalg.norm(a)
    b = np.cross(v, a)
    t = np.linspace(0, 2*np.pi, 721); r = np.deg2rad(angle)
    return np.cos(r)*v[:, None] + np.sin(r)*(a[:, None]*np.cos(t)+b[:, None]*np.sin(t))


class SkyModel:
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

    def _body(self, unix, name):
        meta = self.manifest[name]
        u = (unix-meta['start_unix'])/meta['step_seconds']; i = int(np.floor(u)); f = u-i
        if i < 0 or i >= meta['count']-1:
            raise ValueError('UTC date must be within 1900–2100')
        data = self._array(meta)
        v = data[i].astype(float)*(1-f)+data[i+1].astype(float)*f
        return v/np.linalg.norm(v)

    def evaluate(self, directions, date, band='VIS', model='zodiacal_dgl',
                 sun_avoidance=True, sun_separation=90., moon_avoidance=True, moon_separation=60.):
        """Return unmasked component fluxes plus a common avoidance validity mask."""
        if model not in MODELS or band not in BANDS:
            raise ValueError('Unknown model or band')
        if not all(np.isfinite(a) and 0 <= a <= 180 for a in (sun_separation, moon_separation)):
            raise ValueError('Separation angles must be finite and within 0–180 degrees')
        v = np.asarray(directions, dtype=float)
        if v.shape[-1] != 3 or not np.all(np.isfinite(v)) or not np.allclose(np.linalg.norm(v, axis=-1), 1):
            raise ValueError('Directions must be finite unit vectors with shape (..., 3)')
        unix, utc = utc_time(date)
        sun, moon = self._body(unix, 'solar'), self._body(unix, 'lunar')
        delta = (sun_longitude(sun)-sun_longitude(self.manifest['solar']['reference_vector'])+np.pi)%(2*np.pi)-np.pi
        def sample(data, direction):
            return hp.get_interp_val(data, np.arccos(np.clip(direction[..., 2], -1, 1)),
                                     np.mod(np.arctan2(direction[..., 1], direction[..., 0]), 2*np.pi))
        zodi = sample(self._array(self.manifest['zodiacal'][band]), rotation(v, -delta)) if model != 'dgl' else np.zeros(v.shape[:-1])
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
                   ('MOONOBS', 'GEOCENTRIC'), ('ROT-DEG', np.rad2deg(r['delta'])), ('PHOTCNT', True)]
        if r['model'] != 'dgl': headers += [('REF-UTC', self.manifest['solar']['reference_utc'])]
        if r['model'] != 'zodiacal':
            meta = self.manifest['dgl']['metadata']
            headers += [('CORRMOD', 'C2022_all'), ('CORRSCL', meta['scale']), ('BETAOP', 'MULTIPLY'),
                        ('ALPHAPP', meta['alpha_approximation']), ('BETAREF', meta['beta_reference']), ('DGLFIXED', True)]
        path = Path(filename)
        if not path.name.endswith('.fits.gz'): raise ValueError('FITS output must end with .fits.gz')
        path.parent.mkdir(parents=True, exist_ok=True)
        hp.write_map(str(path), arrays, column_names=names, column_units=units, dtype=np.float32,
                     coord='G' if r.get('coords')=='galactic' else 'C', overwrite=True, extra_header=headers)

    def plot(self, result, target=None, contour_step=.2):
        """Plot the masked sky and Sun/Moon/target positions; return a Figure."""
        import matplotlib.pyplot as plt
        if not any(np.isclose(contour_step, i/10) for i in range(1, 11)):
            raise ValueError('Contour interval must be 0.1, 0.2, ... 1.0 mag')
        r = result; mag = np.where(r['valid'], r['mag'], np.nan)
        finite = mag[np.isfinite(mag)]
        lo = np.floor(finite.min()*2)/2 if len(finite) else 18.
        hi = max(lo+.5, np.ceil(finite.max()*2)/2) if len(finite) else 24.
        fig = plt.figure(figsize=(12, 7))
        wavelengths = self.manifest['bands'][r['band']]
        label={'zodiacal':'Zodiacal light','dgl':'DGL','zodiacal_dgl':'Zodiacal light + DGL'}[r['model']]
        title = f"{r.get('coords','equatorial').title()} | {label} | {r['band']} {wavelengths[0]}–{wavelengths[1]} nm | UTC {r['utc'].replace('T', ' ').rstrip('Z')}"
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore',message='Ignoring specified arguments in this call because figure.*',category=UserWarning)
            projected = hp.mollview(mag, fig=fig.number, coord='G' if r.get('coords')=='galactic' else 'C', cmap='viridis_r', min=lo, max=hi,
                                    badcolor='white', bgcolor='white', title=title, unit='AB mag / arcsec²', return_projected_map=True)
        hp.graticule(dpar=30, dmer=30)
        ax = fig.axes[0]; xmin, xmax, ymin, ymax = ax.proj.get_extent()
        levels = np.round(np.arange(np.floor((lo+1e-9)/contour_step)+1, np.ceil((hi-1e-9)/contour_step))*contour_step, 1)
        if len(levels) and finite.size:
            contours = ax.contour(np.linspace(xmin, xmax, projected.shape[1]), np.linspace(ymin, ymax, projected.shape[0]),
                                  projected, levels=levels, colors='black', linewidths=np.linspace(.6, 2.2, len(levels)), alpha=.7)
            ax.clabel(contours, inline=True, fmt='%.1f', fontsize=8)
        for body, name, marker, color, enabled, angle in [(r['sun'], 'Sun', '*', 'red', r['sun_avoidance'], r['sun_separation']),
                                                         (r['moon'], 'Moon', 'o', '#b56d00', r['moon_avoidance'], r['moon_separation'])]:
            ra, dec = coordinates(from_equatorial(body, r.get('coords','equatorial')))
            if enabled and 0 < angle < 180:
                ring = from_equatorial(boundary(body, angle).T,r.get('coords','equatorial')).T
                hp.projplot(np.rad2deg(np.arctan2(ring[1], ring[0])), np.rad2deg(np.arcsin(np.clip(ring[2], -1, 1))),
                            lonlat=True, color=color, linestyle='--', linewidth=1.4)
            hp.projplot(ra, dec, lonlat=True, marker=marker, color=color, markersize=10, markeredgecolor='black', zorder=12)
            hp.projtext(ra, dec, name, lonlat=True, color=color, fontsize=9, fontweight='bold', ha='right' if name=='Sun' else 'left', va='bottom')
        if target is not None:
            ra, dec = target
            plot_ra, plot_dec = coordinates(from_equatorial(vector(ra,dec), r.get('coords','equatorial')))
            hp.projplot(plot_ra, plot_dec, lonlat=True, marker='x', color='#009cff', markersize=10, markeredgewidth=2, zorder=13)
            hp.projtext(plot_ra, plot_dec, ' Target', lonlat=True, color='#006bbb', fontsize=9)
            p = self.point(ra, dec, r['utc'], band=r['band'], model=r['model'], sun_avoidance=r['sun_avoidance'],
                           sun_separation=r['sun_separation'], moon_avoidance=r['moon_avoidance'], moon_separation=r['moon_separation'])
            label = f"{p['ab_mag_arcsec2']:.3f} mag / arcsec² at (RA, DEC) = ({ra % 360:.2f}°, {dec:.2f}°)"
            if not p['visible']: label += ' · excluded by avoidance (unmasked brightness)'
            fig.text(.5, .02, label, ha='center', fontsize=12, fontweight='bold')
        return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ra', type=float, required=True, help='Equatorial RA in degrees')
    parser.add_argument('--dec', type=float, required=True, help='Equatorial DEC in degrees')
    parser.add_argument('--date', '--utc', required=True, help='UTC ISO date/time, e.g. 2026-10-13T12:00:00Z')
    parser.add_argument('--band', choices=BANDS, default='VIS')
    parser.add_argument('--model', choices=MODELS, default='zodiacal_dgl')
    parser.add_argument('--coords', choices=COORD_SYSTEMS, default='equatorial')
    parser.add_argument('--nside', type=int, default=128)
    parser.add_argument('--data-dir', type=Path, default=DEFAULT_DATA)
    parser.add_argument('--sun-separation', type=float, default=90.)
    parser.add_argument('--moon-separation', type=float, default=60.)
    parser.add_argument('--no-sun-avoidance', action='store_true')
    parser.add_argument('--no-moon-avoidance', action='store_true')
    parser.add_argument('--contour-step', type=float, choices=[i/10 for i in range(1, 11)], default=.2)
    parser.add_argument('--show', action='store_true', help='Open an interactive map window')
    parser.add_argument('--pdf', type=Path, help='Save displayed map as PDF')
    parser.add_argument('--png', type=Path, help='Save displayed map as PNG')
    parser.add_argument('--fits', type=Path, help='Save masked HEALPix map as .fits.gz')
    parser.add_argument('--json', type=Path, help='Save sightline query results as JSON')
    args = parser.parse_args()
    try:
        sky = SkyModel(args.data_dir)
        settings = dict(band=args.band, model=args.model, sun_avoidance=not args.no_sun_avoidance,
                        sun_separation=args.sun_separation, moon_avoidance=not args.no_moon_avoidance,
                        moon_separation=args.moon_separation)
        point = sky.point(args.ra, args.dec, args.date, **settings)
        print(f"{point['ab_mag_arcsec2']:.6f} AB mag / arcsec² at (RA, DEC) = ({point['ra_deg']:.6f}°, {point['dec_deg']:.6f}°)")
        print(f"UTC {point['utc']} | {args.model} | {args.band}")
        print(f"Sun separation {point['sun_separation_deg']:.3f}°; Moon separation {point['moon_separation_deg']:.3f}°; visible: {point['visible']}")
        if not point['visible']: print('Reported brightness is unmasked; this sightline is excluded from the displayed map.')
        if args.json:
            args.json.parent.mkdir(parents=True, exist_ok=True)
            # Strict JSON: non-finite brightness is represented by null.
            clean = {k: (None if isinstance(v, float) and not np.isfinite(v) else v) for k, v in point.items()}
            args.json.write_text(json.dumps(clean, indent=2, allow_nan=False)+'\n')
        if args.show or args.pdf or args.png or args.fits:
            result = sky.sky_map(args.date, args.nside, coords=args.coords, **settings)
            if args.fits: sky.save_fits(result, args.fits)
            if args.show or args.pdf or args.png:
                if not args.show:
                    import matplotlib
                    matplotlib.use('Agg')
                figure = sky.plot(result, (args.ra, args.dec), args.contour_step)
                for output in (args.pdf, args.png):
                    if output:
                        output.parent.mkdir(parents=True, exist_ok=True)
                        figure.savefig(output, dpi=180, bbox_inches='tight')
                import matplotlib.pyplot as plt
                if args.show: plt.show()
                else: plt.close(figure)
    except (ValueError, OSError, KeyError) as exc:
        parser.error(str(exc))


if __name__ == '__main__':
    main()

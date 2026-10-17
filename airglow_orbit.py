#!/usr/bin/env python3
"""Ideal dawn-dusk / noon-midnight airglow maps using bundled browser data.

GLOW-derived calculations are governed by docs/data/glowlicense.txt.
No TLE, orbit perturbations, UV absorption or twilight chemistry history.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
from pathlib import Path
import warnings
import numpy as np
import healpy as hp
from sdrift_sky import SkyModel, DEFAULT_DATA, BANDS, utc_time, coordinates, vector, boundary, ab_magnitude

R = 6371.
ORBIT_TYPES = ('dawn-dusk', 'noon-midnight')
try:
    from numba import njit, prange
    HAVE_NUMBA = True
except ImportError:
    HAVE_NUMBA = False
    prange = range


def orbit_state(sun, height=650., orbit='dawn-dusk', orbital_phase_deg=0., unix=0., epoch=None):
    """North-referenced orbital phase, optionally advanced from an explicit UTC epoch.

    Both orbits: phi=0 projected north, 180 south. Positive direction:
    DD phi=90: s cross north; NM phi=90: midnight, 270 noon.
    Solar-following planes are ideal.
    """
    if orbit not in ORBIT_TYPES or not np.isfinite(height) or not 0 < height < 1500 or not np.isfinite(orbital_phase_deg):
        raise ValueError('Choose an orbit; 0 < H < 1500 km; finite orbital_phase_deg')
    s = np.asarray(sun, dtype=float); s = s/np.linalg.norm(s)
    k = np.array([0., 0., 1.]) if abs(s[2]) < .999999 else np.array([1., 0., 0.])
    north = k-(k@s)*s; north /= np.linalg.norm(north)
    east = np.cross(s, north)
    period = 2*np.pi*np.sqrt((R+height)**3/398600.4418)
    phi = (orbital_phase_deg+(0 if epoch is None else (unix-epoch)*360/period)) % 360
    a = np.deg2rad(phi)
    zenith = np.cos(a)*north+np.sin(a)*(east if orbit == 'dawn-dusk' else -s)
    zenith /= np.linalg.norm(zenith)
    c = np.clip(zenith@s, -1., 1.)
    if abs(abs(c)-1.) < 1e-14: c = np.sign(c)
    z = zenith-c*s
    z = z/np.linalg.norm(z) if np.linalg.norm(z) > 1e-10 else north
    # Columns map the reference Sun +X / satellite in XZ to the real sky.
    axes = np.column_stack((s, np.cross(z, s), z))
    return dict(sun=s, zenith=zenith, sza=float(np.rad2deg(np.arccos(c))),
                orbital_phase_deg=float(phi), period=float(period), axes=axes,
                limb_zenith_deg=float(180-np.rad2deg(np.arcsin(R/(R+height)))))


def _kernel(height, sza, directions, q, z, chi, support, step):
    """Midpoint LOS integral, same samples/units as the browser integrator."""
    out = np.full((len(directions), 6), np.nan)
    rr = R+height; angle = np.deg2rad(sza)
    px, pz = rr*np.cos(angle), rr*np.sin(angle)
    for i in prange(len(directions)):
        v = directions[i]; b = px*v[0]+pz*v[2]
        if b < 0 and b*b-rr*rr+R*R >= 0:
            out[i] = (np.nan, np.nan, 1., 0., 0., 0.)
            continue
        hmin = np.sqrt(max(rr*rr-min(b, 0.)**2, 0.))-R
        disc = b*b-rr*rr+(R+z[-1])**2
        root = np.sqrt(max(0., disc)); start = max(-b-root, 0.); end = max(-b+root, 0.) if disc >= 0 else 0.
        count = max(1, int(np.ceil((end-start)/step))); ds = (end-start)/count
        flux = 0.; lit = 0.; cutoff = False
        for j in range(count):
            t = start+ds*(j+.5); x = px+t*v[0]; y = t*v[1]; zz = pz+t*v[2]
            rad = np.sqrt(x*x+y*y+zz*zz); h = rad-R
            if h < z[0] or h > z[-1]: continue
            a = np.arccos(min(1., max(-1., x/rad)))*180/np.pi
            shadow = x < 0 and y*y+zz*zz <= R*R
            if not shadow and a >= 1.85*180/np.pi: cutoff = True
            iz = min(len(z)-2, max(0, np.searchsorted(z, h, side='right')-1))
            ic = min(len(chi)-2, max(0, np.searchsorted(chi, a, side='right')-1))
            fz = (h-z[iz])/(z[iz+1]-z[iz]); fc = (a-chi[ic])/(chi[ic+1]-chi[ic])
            value = (1-fc)*((1-fz)*q[ic,iz,0]+fz*q[ic,iz+1,0])+fc*((1-fz)*q[ic+1,iz,0]+fz*q[ic+1,iz+1,0])
            if shadow:
                value += (1-fc)*((1-fz)*q[ic,iz,1]+fz*q[ic,iz+1,1])+fc*((1-fz)*q[ic+1,iz,1]+fz*q[ic+1,iz+1,1])
            flux += value
            if not shadow: lit += value
        out[i] = (flux*ds*1e5/(4*np.pi), lit*ds*1e5/(4*np.pi), 0., (1. if hmin >= support and not cutoff else 0.), hmin, (1. if cutoff else 0.))
    return out


if HAVE_NUMBA:
    _fast_kernel = njit(parallel=True, cache=True)(_kernel)


def _numpy_kernel(height, sza, directions, q, z, chi, support, step):
    """NumPy fallback: vectorize all midpoint samples within each sightline."""
    from math import cos, sin, radians
    pos = (R+height)*np.array([cos(radians(sza)), 0., sin(radians(sza))])
    result = np.full((len(directions), 6), np.nan)
    for i, v in enumerate(directions):
        b = v@pos; rr = (R+height)**2
        if b < 0 and b*b-rr+R*R >= 0:
            result[i] = [np.nan, np.nan, 1, 0, 0, 0]; continue
        hmin = np.sqrt(max(rr-min(b, 0)**2, 0))-R
        disc = b*b-rr+(R+z[-1])**2; root = np.sqrt(max(0, disc))
        start = max(-b-root, 0); end = max(-b+root, 0) if disc >= 0 else 0
        count = max(1, int(np.ceil((end-start)/step))); ds = (end-start)/count
        points = pos+(start+ds*(np.arange(count)+.5))[:, None]*v
        radii = np.linalg.norm(points, axis=1); h = radii-R
        keep = (h >= z[0]) & (h <= z[-1]); points = points[keep]; h = h[keep]; radii = radii[keep]
        a = np.rad2deg(np.arccos(np.clip(points[:, 0]/radii, -1, 1)))
        shadow = (points[:, 0] < 0) & (np.sum(points[:, 1:]**2, axis=1) <= R*R)
        cutoff = np.any((~shadow) & (a >= np.rad2deg(1.85)))
        iz = np.clip(np.searchsorted(z, h, side='right')-1, 0, len(z)-2)
        ic = np.clip(np.searchsorted(chi, a, side='right')-1, 0, len(chi)-2)
        fz = ((h-z[iz])/(z[iz+1]-z[iz]))[:, None]; fc = ((a-chi[ic])/(chi[ic+1]-chi[ic]))[:, None]
        values = (1-fc)*((1-fz)*q[ic,iz]+fz*q[ic,iz+1])+fc*((1-fz)*q[ic+1,iz]+fz*q[ic+1,iz+1])
        emission = values[:, 0]+np.where(shadow, values[:, 1], 0.)
        factor = ds*1e5/(4*np.pi)
        result[i] = [emission.sum()*factor, emission[~shadow].sum()*factor, 0, hmin >= support and not cutoff, hmin, cutoff]
    return result


class AirglowOrbitModel(SkyModel):
    def evaluate(self, directions, date, band='VIS', model='combined', height=650., orbit='dawn-dusk', orbital_phase_deg=0.,
                 epoch=None, sun_avoidance=True, sun_separation=90.,
                 moon_avoidance=True, moon_separation=60., engine='auto', step_km=2., phase=None, min_elv_deg=0., kind='airglow'):
        if kind not in ('airglow', 'total'): raise ValueError('Kind must be airglow or total')
        if phase is not None: orbital_phase_deg = phase  # alias; NORTH_V1 convention
        if band not in BANDS or model not in ('glow', 'combined'): raise ValueError('Unknown band/model')
        if not all(np.isfinite(x) for x in (sun_separation, moon_separation, step_km)) or not (0 <= sun_separation <= 180 and 0 <= moon_separation <= 180 and step_km > 0):
            raise ValueError('Separation angles must be 0-180; positive LOS step')
        if not np.isfinite(min_elv_deg) or not -90 <= min_elv_deg <= 180:
            raise ValueError('Minimum ELV must be finite and between -90 and 180 degrees')
        if engine not in ('auto', 'numba', 'numpy'): raise ValueError('Unknown engine')
        if engine == 'numba' and not HAVE_NUMBA: raise ValueError('Install numba or use --engine numpy')
        v = np.asarray(directions, dtype=float).reshape(-1, 3)
        if not np.isfinite(v).all() or not np.allclose(np.linalg.norm(v, axis=1), 1.): raise ValueError('LOS must be finite unit vectors')
        unix, utc = utc_time(date); epoch_unix = None if epoch is None else utc_time(epoch)[0]
        sun, moon = self._body(unix, 'solar'), self._body(unix, 'lunar')
        state = orbit_state(sun, height, orbit, orbital_phase_deg, unix, epoch_unix)
        theta = np.rad2deg(np.arccos(np.clip(v@state['zenith'], -1, 1)))
        elv = state['limb_zenith_deg']-theta
        view = elv >= min_elv_deg-1e-9
        if sun_avoidance: view &= v@sun <= np.cos(np.deg2rad(sun_separation))+1e-12
        if moon_avoidance: view &= v@moon <= np.cos(np.deg2rad(moon_separation))+1e-12
        info = self.manifest['airglow'][model]; q = self._array(info['bands'][band]).astype(float)
        function = _fast_kernel if HAVE_NUMBA and engine != 'numpy' else _numpy_kernel
        a = np.full((len(v), 6), np.nan); a[:, [2, 3, 5]] = 0
        a[view] = function(height, state['sza'], np.ascontiguousarray(v[view]@state['axes']), q,
                           np.array(info['altitude_km']), np.array(info['chi_deg']),
                           float(info['metadata'].get('support_min_altitude_km', info['altitude_km'][0])), step_km)
        lo, hi = self.manifest['bands'][band]
        with np.errstate(divide='ignore', invalid='ignore'):
            mag = -2.5*np.log10(6.62607015e-27*a[:, 0]/(np.log(hi/lo)*206264.80624709636**2)/(3631e-23))
        result = dict(flux=a[:, 0], lit=a[:, 1], earth=a[:, 2].astype(bool), valid=a[:, 3].astype(bool),
                    hmin=a[:, 4], cutoff=a[:, 5].astype(bool), theta=theta, elv=elv, view=view, mag=mag, min_elv_deg=min_elv_deg,
                    state=state, sun=sun, moon=moon, utc=utc, band=band, model=model, height=height,
                    orbit=orbit, orbital_phase_deg=orbital_phase_deg, epoch=epoch,
                    sun_avoidance=sun_avoidance, sun_separation=sun_separation,
                    moon_avoidance=moon_avoidance, moon_separation=moon_separation, step_km=step_km,
                    engine='numba' if function is globals().get('_fast_kernel') else 'numpy', kind=kind)
        if kind == 'total':
            background = SkyModel.evaluate(self, v, date, band=band, model='zodiacal_dgl',
                sun_avoidance=sun_avoidance, sun_separation=sun_separation,
                moon_avoidance=moon_avoidance, moon_separation=moon_separation)
            shown = view & ~result['earth']
            result['airglow_photon_radiance'] = result['flux'].copy()
            result['airglow'] = 6.62607015e-10*result['flux']/np.log(hi/lo)  # MJy/sr
            for key in ['zodiacal', 'dgl']:
                result[key] = np.where(shown, background[key], np.nan)
            result['flux'] = result['zodiacal']+result['dgl']+result['airglow']
            result['mag'] = ab_magnitude(result['flux'])
            result['delta'] = background['delta']
        return result

    def point(self, ra, dec, date, **settings):
        r = self.evaluate(vector(ra, dec), date, **settings)
        point = dict(utc=r['utc'], ra_deg=ra % 360, dec_deg=dec, band=r['band'], orbit=r['orbit'],
                    height_km=r['height'], orbital_phase_deg=r['state']['orbital_phase_deg'], sza_deg=r['state']['sza'],
                    zenith_deg=float(r['theta'][0]), elv_deg=float(r['elv'][0]), minimum_elv_deg=r['min_elv_deg'], view_allowed=bool(r['view'][0]),
                    earth_occulted=bool(r['earth'][0]), model_valid=bool(r['valid'][0]),
                    ab_mag_arcsec2=float(r['mag'][0]), photon_radiance=float(r.get('airglow_photon_radiance',r['flux'])[0]), kind=r['kind'],
                    sun_ra_dec_deg=coordinates(r['sun']), moon_ra_dec_deg=coordinates(r['moon']),
                    sun_separation_deg=float(np.rad2deg(np.arccos(np.clip(vector(ra, dec)@r['sun'], -1, 1)))),
                    moon_separation_deg=float(np.rad2deg(np.arccos(np.clip(vector(ra, dec)@r['moon'], -1, 1)))))

        if r['kind'] == 'total':
            point.update(intensity_mjy_sr=float(r['flux'][0]), zodiacal_mjy_sr=float(r['zodiacal'][0]),
                         dgl_mjy_sr=float(r['dgl'][0]), airglow_mjy_sr=float(r['airglow'][0]))
        return point

    def save_fits(self, r, filename):
        path = Path(filename)
        if not path.name.endswith('.fits.gz'): raise ValueError('FITS output must end in .fits.gz')
        path.parent.mkdir(parents=True, exist_ok=True)
        keys = ['mag', 'flux', 'lit', 'earth', 'valid', 'hmin', 'cutoff', 'theta', 'view', 'elv']
        names = ['AB_MAG_ARCSEC2', 'PHOTON_RADIANCE', 'SUNLIT_RADIANCE', 'EARTH_OCCULTED', 'MODEL_VALID', 'MIN_ALTITUDE', 'PHOTO_CUTOFF', 'ZENITH_ANGLE', 'VIEW_ALLOWED', 'ELV']
        units = ['mag/arcsec2', 'photons/s/cm2/sr', 'photons/s/cm2/sr', '', '', 'km', '', 'deg', '', 'deg']
        if r['kind'] == 'total':
            names[1] = 'INTENSITY_MJY_SR'; units[1] = 'MJy/sr'
            keys += ['zodiacal', 'dgl', 'airglow', 'airglow_photon_radiance']
            names += ['ZODIACAL_MJY_SR', 'DGL_MJY_SR', 'AIRGLOW_MJY_SR', 'AIRGLOW_PHOTON_RADIANCE']
            units += ['MJy/sr', 'MJy/sr', 'MJy/sr', 'photons/s/cm2/sr']
        lo, hi = self.manifest['bands'][r['band']]
        headers = [('DATE-UTC', r['utc']), ('ORBIT', r['orbit']), ('ORBMODEL', 'IDEAL SOLAR CIRCLE'),
                   ('SAT_H', r['height']), ('SAT_SZA', r['state']['sza']), ('PHASE', r['state']['orbital_phase_deg']),
                   ('PHREF', 'NORTH'), ('PHCONV', 'NORTH_V1'), ('PHASE0', r['orbital_phase_deg']), ('PHMODE', 'EPOCH' if r['epoch'] else 'MANUAL'),
                   ('PERIOD', r['state']['period']), ('ELVMIN', r['min_elv_deg']), ('ELVREF', 'SPHERICAL EARTH'), ('EARTHR', R),
                   ('SUNMASK', r['sun_avoidance']), ('ELONGMIN', r['sun_separation']),
                   ('MOONMASK', r['moon_avoidance']), ('MOONMIN', r['moon_separation']), ('MOONOBS', 'GEOCENTRIC'),
                   ('BAND', r['band']), ('BANDMIN', lo/1000), ('BANDMAX', hi/1000), ('PHOTCNT', True),
                   ('FULLBAND', False), ('ABSORPT', False), ('TWILHIST', False), ('INTVER', 2), ('DS_KM', r['step_km']),
                   ('TABLESHA', self.manifest['airglow'][r['model']]['source_sha256']), ('MODEL', r['model'])]
        headers.append(('DATATYPE', r['kind'].upper()))
        if r['kind'] == 'total':
            meta = self.manifest['dgl']['metadata']
            headers += [('REF-UTC', self.manifest['solar']['reference_utc']), ('ROT-DEG', np.rad2deg(r['delta'])),
                        ('SRC_SHA', self.manifest['zodiacal'][r['band']]['source_sha256']),
                        ('DGL_SHA', self.manifest['dgl']['bands'][r['band']]['source_sha256']),
                        ('CORRMOD','C2022_all'), ('CORRSCL',meta['scale']), ('BETAOP','MULTIPLY'),
                        ('ALPHAPP',meta['alpha_approximation']), ('BETAREF',meta['beta_reference']), ('DGLFIXED',True)]
        if r['epoch']: headers.append(('EPOCH', r['epoch']))
        for body, prefix in [(r['sun'], 'SUN'), (r['moon'], 'MOON'), (r['state']['zenith'], 'ZEN')]:
            ra, dec = coordinates(body); headers += [(prefix+'RA', ra), (prefix+'DEC', dec)]
        hp.write_map(str(path), [r[k] for k in keys], column_names=names, column_units=units, dtype=np.float32,
                     coord='C', overwrite=True, extra_header=headers)

    def plot(self, r, contour_step=.5, figure=None):
        import matplotlib.pyplot as plt
        with plt.rc_context({'text.usetex': False, 'font.family': 'DejaVu Sans'}):
            return self._plot(r, contour_step, figure)

    def _plot(self, r, contour_step=.5, figure=None):
        import matplotlib.pyplot as plt
        if not any(np.isclose(contour_step, i/10) for i in range(1, 11)): raise ValueError('Contour step must be 0.1 ... 1.0')
        fig = figure if figure is not None else plt.figure(figsize=(13, 7.5))
        fig.clear(); ax = fig.add_subplot(111, projection='mollweide')
        x = np.linspace(-np.pi, np.pi, 721); y = np.linspace(-np.pi/2, np.pi/2, 361)
        xx, yy = np.meshgrid(x, y)
        ids = hp.ang2pix(r['nside'], np.pi/2-yy, np.mod(-xx, 2*np.pi))
        mag = r['mag'][ids]; finite = r['mag'][np.isfinite(r['mag'])]
        lo = np.floor(finite.min()*2)/2 if len(finite) else 18.
        hi = max(lo+.5, np.ceil(finite.max()*2)/2) if len(finite) else 24.
        cmap = plt.get_cmap('viridis_r').copy(); cmap.set_bad('#f5f5f5')
        im = ax.pcolormesh(x, y, np.ma.masked_invalid(mag), shading='auto', cmap=cmap, vmin=lo, vmax=hi, rasterized=True)
        earth = r['earth'][ids] & r['view'][ids]
        ax.contourf(xx, yy, earth, levels=[.5, 1.5], colors=['#aeb5be'])
        unsupported = r['view'][ids] & ~r['earth'][ids] & ~r['valid'][ids]
        ax.contourf(xx, yy, unsupported, levels=[.5, 1.5], colors='none', hatches=['///'])
        levels = np.round(np.arange(np.floor((lo+1e-9)/contour_step)+1, np.ceil((hi-1e-9)/contour_step))*contour_step, 1)
        if len(finite) and len(levels):
            good = np.ma.masked_where(~r['valid'][ids] | ~np.isfinite(mag), mag)
            if good.count():
                cs = ax.contour(xx, yy, good, levels=levels, colors='black', linewidths=np.linspace(.6, 2.2, len(levels)), alpha=.65)
                ax.clabel(cs, fmt='%.1f', fontsize=8)
        def project(v):
            return -(np.arctan2(v[1], v[0])+np.pi)%(2*np.pi)-np.pi, np.arcsin(np.clip(v[2], -1, 1))
        def ring(body, angle, color):
            if not 0 < angle < 180: return
            lon, lat = project(boundary(body, angle)); lon = lon.copy()
            lon[1:][np.abs(np.diff(lon)) > np.pi] = np.nan
            ax.plot(lon, lat, '--', color=color, lw=1.3)
        for v, name, marker, color, enabled, angle in [(r['sun'], 'Sun', '*', 'red', r['sun_avoidance'], r['sun_separation']),
                  (r['moon'], 'Moon', 'o', '#b56d00', r['moon_avoidance'], r['moon_separation'])]:
            if enabled: ring(v, angle, color)
            a, b = project(v); ax.plot(a, b, marker, color=color, ms=9, markeredgecolor='black'); ax.annotate(name, (a,b), xytext=(6,6), textcoords='offset points', color=color, weight='bold', fontsize=9)
        ring(r['state']['zenith'], r['state']['limb_zenith_deg']-r['min_elv_deg'], '#008060')
        a,b = project(r['state']['zenith']); ax.plot(a,b,'o',color='#1579bf',ms=4); ax.annotate('Zenith',(a,b),xytext=(6,-12),textcoords='offset points',color='#1579bf',fontsize=9)
        ax.grid(alpha=.3); ticks = np.arange(-150, 180, 30)
        ax.set_xticks(np.deg2rad(ticks), [f'{int(-t)%360}°' for t in ticks]); ax.set_xlabel('RA (deg; increases to the left)'); ax.set_ylabel('DEC (deg)')
        lo_nm, hi_nm = self.manifest['bands'][r['band']]
        label = 'Zodiacal light + DGL + Airglow' if r['kind'] == 'total' else 'Airglow'
        fig.suptitle(f"{label} | {r['band']} {lo_nm:g}-{hi_nm:g} nm | UTC {r['utc'].replace('T', ' ').rstrip('Z')}\n"
                     f"{r['orbit']} | H={r['height']:g} km | φ={r['state']['orbital_phase_deg']:.1f}° | SZA={r['state']['sza']:.1f}° | ELV ≥ {r['min_elv_deg']:g}°", fontsize=11)
        fig.colorbar(im, ax=ax, orientation='horizontal', fraction=.045, pad=.12, label='AB mag / arcsec²')
        fig.text(.5,.025, f"Sun separation ≥ {r['sun_separation']:g}° ({'on' if r['sun_avoidance'] else 'off'}); Moon ≥ {r['moon_separation']:g}° ({'on' if r['moon_avoidance'] else 'off'}) | ELV ≥ {r['min_elv_deg']:g}° | Gray: Earth; hatch: partial model; pale: excluded\nIdeal solar-following circle; mean emission table; no UV absorption or twilight history", ha='center', fontsize=8)
        return fig


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gui', action='store_true', help='Interactive Python controls')
    p.add_argument('--date', '--utc', default=datetime.now(timezone.utc).isoformat())
    p.add_argument('--orbit', choices=ORBIT_TYPES, default='dawn-dusk')
    p.add_argument('--height', '--H', type=float, default=650.)
    p.add_argument('--orbital-phase-deg', '--phase', dest='orbital_phase_deg', type=float, default=0.,
                   help='North-referenced orbital phase phi, in degrees; 0=north, 180=south')
    p.add_argument('--epoch', help='UTC of the input north phase; omit to hold phase fixed')
    p.add_argument('--min-elv', '--elv-min', type=float, default=0., help='Minimum Earth limb elevation, deg (default: 0)')
    p.add_argument('--band', choices=BANDS, default='VIS')
    p.add_argument('--kind', choices=['airglow','total'], default='airglow', help='Total = zodiacal + DGL + airglow')
    p.add_argument('--model', choices=['combined','glow'], default='combined')
    p.add_argument('--nside', type=int, default=128)
    p.add_argument('--data-dir', type=Path, default=DEFAULT_DATA)
    p.add_argument('--engine', choices=['auto','numba','numpy'], default='auto')
    p.add_argument('--sun-separation', type=float, default=90.)
    p.add_argument('--moon-separation', type=float, default=60.)
    p.add_argument('--no-sun-avoidance', action='store_true')
    p.add_argument('--no-moon-avoidance', action='store_true')
    p.add_argument('--contour-step', type=float, choices=[i/10 for i in range(1,11)], default=.5)
    p.add_argument('--ra', type=float); p.add_argument('--dec', type=float)
    p.add_argument('--show', action='store_true')
    for kind in ['pdf','png','fits']: p.add_argument('--'+kind, type=Path)
    a = p.parse_args()
    if a.gui:
        from airglow_orbit_gui import main as gui
        gui(a.data_dir, kind=a.kind); return
    try:
        sky = AirglowOrbitModel(a.data_dir)
        settings = dict(band=a.band, model=a.model, kind=a.kind, height=a.height, orbit=a.orbit, orbital_phase_deg=a.orbital_phase_deg, epoch=a.epoch,
                        min_elv_deg=a.min_elv, sun_avoidance=not a.no_sun_avoidance,
                        sun_separation=a.sun_separation, moon_avoidance=not a.no_moon_avoidance,
                        moon_separation=a.moon_separation, engine=a.engine)
        if (a.ra is None) != (a.dec is None): p.error('Provide both --ra and --dec')
        if a.ra is not None: print(sky.point(a.ra,a.dec,a.date,**settings))
        if a.pdf or a.png or a.fits or a.show or a.ra is None:
            r = sky.sky_map(a.date, nside=a.nside, **settings)
            print(f"{r['orbit']} | UTC {r['utc']} | north phase {r['state']['orbital_phase_deg']:.3f} deg | SZA {r['state']['sza']:.3f} deg | {r['engine']} | {int(r['view'].sum())} allowed pixels")
            if a.fits: sky.save_fits(r,a.fits)
            if a.pdf or a.png or a.show:
                if not a.show:
                    import matplotlib
                    matplotlib.use('Agg')
                fig = sky.plot(r,contour_step=a.contour_step)
                for path in [a.pdf,a.png]:
                    if path: path.parent.mkdir(parents=True,exist_ok=True); fig.savefig(path,dpi=160)
                if a.show:
                    import matplotlib.pyplot as plt
                    plt.show()
    except (ValueError,OSError) as exc: p.exit(2,str(exc)+'\n')

if __name__ == '__main__': main()

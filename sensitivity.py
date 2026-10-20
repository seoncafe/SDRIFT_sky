"""Gaussian surface-brightness detection limit, using S-DRIFT sky backgrounds.

Notebook CCD equation plus a fully correlated sky-only systematic residual.
Exposure is usable accumulated time at the specified (fixed) sky geometry.
No Earthshine, HST correction, filter leakage or orbit duty-cycle integration.
"""
from __future__ import annotations
import numpy as np

BAND_EDGES_NM = {'UV1': (200.,400.), 'UV2': (250.,400.),
                 'UV3': (300.,400.), 'VIS': (400.,700.)}
DEFAULTS = dict(f_sys=.001, exposure_hours=10., bin_side=5, snr=3.,
                diameter_cm=30., pixel_scale_arcsec=1.96, read_noise_e=1.,
                dark_current_e_s=.05, n_reads=1, throughput=None)

def parameters(band='VIS', **options):
    """Validate instrument settings; f_sys is a fraction, not a percentage."""
    if band not in BAND_EDGES_NM: raise ValueError('Unknown sensitivity band')
    unknown = set(options)-set(DEFAULTS)
    if unknown: raise TypeError(f'Unknown sensitivity settings: {sorted(unknown)}')
    p = DEFAULTS | options
    if p['throughput'] is None: p['throughput'] = .9 if band == 'VIS' else .5
    for k,v in p.items():
        if isinstance(v, (bool, np.bool_)) or not np.isscalar(v) or not np.isfinite(v):
            raise ValueError(f'{k} must be finite and numeric')
        if k in ('f_sys','read_noise_e','dark_current_e_s'):
            if v < 0: raise ValueError(f'{k} must be nonnegative')
        elif v <= 0: raise ValueError(f'{k} must be positive')
    if p['f_sys'] > 1 or p['throughput'] > 1:
        raise ValueError('f_sys and throughput must be fractions within 0–1')
    for k in ('bin_side','n_reads'):
        if int(p[k]) != p[k]: raise ValueError(f'{k} must be an integer')
        p[k] = int(p[k])
    return p

def calculate_sensitivity(background_mjy_sr, band='VIS', **options):
    """Return limiting AB mag/arcsec² and count budget for scalar/array sky.

    N = bin_side**2 detector pixels (independent of HEALPix NSIDE), binned in
    software after readout: each pixel adds RN**2 per read.  On-chip binning
    would add RN**2 once per read (negligible for sky-limited bins).
    S = SNR * sqrt(S + B_sky + B_dark + N*n_reads*RN**2 + (f_sys*B_sky)**2).
    A band-mean photon-counting MJy/sr brightness exactly determines the counts
    for flat throughput. The systematic is correlated across time and bin.
    NaN background remains NaN, e.g. masked or unsupported sightlines.
    """
    p = parameters(band, **options)
    bg = np.asarray(background_mjy_sr, dtype=float)
    if np.any(np.isinf(bg)) or np.any(bg < 0):
        raise ValueError('Background must be nonnegative and finite, or NaN for masked pixels')
    lo,hi = BAND_EDGES_NM[band]
    area = np.pi*(p['diameter_cm']/2)**2
    # 1 MJy/sr = 1e-17 erg s^-1 cm^-2 Hz^-1 sr^-1; integrate dnu/nu.
    rate_per_mjy = (1e-17 / 6.62607015e-27 * np.log(hi/lo) * area
                    * (p['pixel_scale_arcsec']/206264.80624709636)**2 * p['throughput'])
    t = p['exposure_hours']*3600
    n = p['bin_side']**2
    sky = bg*rate_per_mjy*n*t
    dark = n*p['dark_current_e_s']*t
    read_var = n*p['n_reads']*p['read_noise_e']**2
    sys = p['f_sys']*sky
    variance = sky+dark+read_var+sys**2
    source = .5*p['snr']**2*(1+np.sqrt(1+4*variance/p['snr']**2))
    limit = source/(rate_per_mjy*n*t)
    with np.errstate(invalid='ignore', divide='ignore'):
        mag = -2.5*np.log10(limit*1e6/206264.80624709636**2/3631.)
    return dict(limiting_ab_mag_arcsec2=mag, limiting_mjy_sr=limit,
                background_mjy_sr=bg, sky_counts_e=sky, dark_counts_e=dark,
                read_variance_e2=read_var, systematic_sigma_e=sys,
                source_counts_e=source, rate_per_mjy_sr_e_s_pixel=rate_per_mjy,
                systematic_floor_mjy_sr=p['snr']*p['f_sys']*bg, sensitivity=p)

# S-DRIFT Python interface and usage

Status: 2026-10-20, common Python API v2.

Zodiacal light, DGL, airglow, and their sum are computed through **one
`SkyModel` interface**. All modes share the same function names, settings,
intensity units, and observing conditions. This document describes the Python
API, the command line (CLI), the GUI, and the saved files.

## 1. Quick start

```python
from sdrift_sky import SkyModel

sky = SkyModel()
point = sky.point(
    0, 60, "2026-10-13T12:00:00Z",
    model="zodiacal_dgl",
    band="VIS",
)
print(point["ab_mag_arcsec2"])
print(point["visible"])
```

In `point(0, 60, date, ...)`, **0 is RA = 0° and 60 is DEC = +60°**. RA is in
degrees, not hours; RA = 0° corresponds to 0 h. The trailing `Z` of the date
means UTC. The example computes the VIS-band zodiacal light + DGL brightness in
that direction at that UTC.

The orbit, ELV, and Sun/Moon avoidance conditions are applied by default. For
directions with `visible=False`, brightness and intensity are `NaN`, so always
check the visibility flag together with the brightness.

To change the mode, **change only `model`**.

```python
from sdrift_sky import SkyModel

sky = SkyModel()
for model in ["zodiacal", "dgl", "zodiacal_dgl", "airglow", "total"]:
    p = sky.point(
        30, 10, "2026-10-17T00:00:00Z",
        model=model,
        orbit="noon-midnight",
        orbital_phase_deg=90,
    )
    print(model, p["ab_mag_arcsec2"], p["intensity_mjy_sr"], p["visible"])
```

## 2. Required files and installation

After downloading the GitHub repository, keep the following layout.

```text
SDRIFT_sky/
├── sdrift_sky.py               # common SkyModel, CLI
├── airglow_orbit.py            # numerical engine, compatible CLI/class name
├── airglow_orbit_gui.py        # common Python GUI
├── sky_coordinates.py          # Equatorial / Galactic rotation
├── requirements-python.txt
└── docs/data/
    ├── manifest.json           # file names, units, array shapes, sources and hashes
    └── ...                     # all compressed data files
```

Do not copy a single Python file on its own; keep the files above together with
`docs/data/`. The large original FITS/NPZ files, the GLOW/PALACE sources, and a
compiler are not needed to use the published data. No internet connection is
needed at run time.

```sh
python3 -m pip install -r requirements-python.txt
# Optional: install a compatible Numba for parallel airglow integration
python3 -m pip install numba
```

Python 3.10 or later is targeted. NumPy, healpy, Astropy, and Matplotlib are
required; Numba is optional. Use a Numba release compatible with the installed
Python and NumPy. The GUI needs a graphical desktop, Tkinter, and the Tcl/Tk
matching that Python.

```sh
python3 -m tkinter
```

Check that this command opens a Tkinter window. Tkinter is installed as a
package of your Python distribution; do not use `pip install tkinter`.
On a headless server, save files instead of using the GUI or `--show`.

To use a different data folder:

```python
from sdrift_sky import SkyModel
sky = SkyModel(data_dir="/absolute/path/to/docs/data")
```

Data are read on first use, checked against their size and SHA-256, and cached
in the instance. Reuse the `sky` instance when computing many directions or
dates.

## 3. Common modes and defaults

### 3.1 `model`: sky components to sum

| Value | Components |
|---|---|
| `"zodiacal"` | Zodiacal light |
| `"dgl"` | Diffuse Galactic Light |
| `"zodiacal_dgl"` | Zodiacal light + DGL |
| `"airglow"` | Airglow |
| `"total"` | Zodiacal light + DGL + Airglow; default |

All modes use the same orbit and observing conditions. The zodiacal/DGL modes
also restrict the observable region, but they do not compute or add airglow.

### 3.2 `airglow_model`: airglow emission set

| Value | Meaning |
|---|---|
| `"combined"` | GLOW + assumed PALACE night layers; default |
| `"glow"` | GLOW selected emission lines |

`model` and `airglow_model` are independent settings.

```python
from sdrift_sky import SkyModel
sky = SkyModel()
p = sky.point(
    30, 10, "2026-10-17T00:00:00Z",
    model="total", airglow_model="glow",
    orbit="noon-midnight", orbital_phase_deg=90,
)
```

With `model="zodiacal"`, `"dgl"`, or `"zodiacal_dgl"`, `airglow_model` is not
used. The GUI disables the airglow selection in these modes.

### 3.3 Common settings

`point()`, `evaluate()`, and `sky_map()` accept the following settings.

| Python argument | Default | Units, allowed values, meaning |
|---|---|---|
| `model` | `"total"` | one of the five sky modes above |
| `airglow_model` | `"combined"` | `combined`, `glow` |
| `band` | `"VIS"` | `UV1`, `UV2`, `UV3`, `VIS` |
| `orbit` | `"dawn-dusk"` | `dawn-dusk`, `noon-midnight` |
| `height` | `650.0` | km; `0 < height < 1500` |
| `orbital_phase_deg` | `0.0` | phase from the northernmost point, degrees; taken modulo 360° |
| `epoch` | `None` | reference UTC of the phase; `None` keeps the phase fixed |
| `min_elv_deg` | `30.0` | degrees; −90 to 180 |
| `sun_avoidance` | `True` | apply Sun avoidance |
| `sun_separation` | `90.0` | degrees; 0 to 180 |
| `moon_avoidance` | `True` | apply Moon avoidance |
| `moon_separation` | `60.0` | degrees; 0 to 180 |
| `engine` | `"auto"` | `auto`, `numba`, `numpy` |
| `step_km` | `2.0` | maximum airglow integration step, km; positive; API only |

`coords` defaults to `"equatorial"`; `"galactic"` is also available. See
Section 5 for coordinate input and display rules. `nside` is an argument of
`sky_map()` only, and `contour_step` of `plot()` only. `data_dir` is given when
constructing `SkyModel()`.

### 3.4 Bands and brightness units

| Band | Wavelength range |
|---|---|
| `UV1` | 200–400 nm = 0.2–0.4 μm |
| `UV2` | 250–400 nm = 0.25–0.4 μm |
| `UV3` | 300–400 nm = 0.3–0.4 μm |
| `VIS` | 400–700 nm = 0.4–0.7 μm |

All bands assume a constant photon-counting throughput.
**Intensities are in MJy/sr in every mode**, and brightness is in
**AB mag/arcsec²**. Airglow-only mode uses the same units. The separate
`photon_radiance` and `airglow_photon_radiance` fields hold the airglow value in
photon units. Components that are not selected are 0 in observable directions.

Sums are formed by adding intensities in the same band and then converting to
AB magnitude; AB magnitudes are never added directly. A larger AB magnitude is
fainter.

## 4. Date, orbit, and observing conditions

### 4.1 UTC

Pass the date as an ISO string.

```text
2026-10-19T00:00:00Z
2026-10-19T09:00:00+09:00
```

These two times are identical. Strings without a time zone are accepted but
interpreted as UTC, so giving `Z` or an offset explicitly is safer. The stored
Sun and Moon data cover 1900–2100. `date` is required in the API; the CLI/GUI
use the current UTC when it is omitted.

### 4.2 Orbit model

Both orbits are circular sun-synchronous orbits of radius R + H (R = 6371 km,
spherical Earth) with J2 secular rates only; there is no TLE/SGP4 propagation,
drag, or higher harmonics.

- The inclination follows from the J2 sun-synchronous condition
  dΩ/dt = −1.5 n J2 (R_eq/a)² cos i = 360° per tropical year
  (J2 = 1.08262668e-3, R_eq = 6378.137 km, μ = 398600.4418 km³/s²).
  At H = 650 km, i = 97.96°.
- The ascending node is at a fixed local time in **mean** solar time:
  LTAN 18:00 for dawn-dusk and 12:00 for noon-midnight.
  RAAN = mean-Sun RA (of date) + 15° × (LTAN − 12 h). The plane therefore
  follows the mean Sun; the true Sun differs by the equation of time and by its
  declination.
- Earth's pole is the mean pole of date (IAU 1976 precession; agrees with ERFA
  to < 0.3 arcsec over 1900–2100). Nutation (< 0.005°) and the frame bias are
  neglected.

As a result the **satellite SZA varies around the orbit**. For dawn-dusk at
650 km it spans about 82–98° at the equinoxes, 74.5–105.5° at the December
solstice, and 58.6–121.4° at the June solstice, so near the June solstice the
satellite enters Earth's shadow (SZA > 114.85° at 650 km). For noon-midnight
the SZA spans about 2–178°.

### 4.3 Orbital phase φ (convention NORTH_V2)

φ is measured along the orbit from its northernmost point and advances in the
direction of motion.

| Phase | Dawn-dusk | Noon-midnight |
|---|---|---|
| 0° | northernmost point | northernmost point |
| 90° | descending node, dawn side | descending node, midnight side |
| 180° | southernmost point | southernmost point |
| 270° | ascending node, dusk side | ascending node, noon side |

"Northernmost" refers to the orbit itself (the point of maximum declination,
about 82° for i = 97.96°), not to a position above the geographic pole.

### 4.4 Fixed phase and time evolution

With `epoch=None`, the given phase is kept when the date changes; the Sun, the
Moon, and the orbit plane still move with the date. If `epoch` is given, the
satellite is assumed to be at the given phase at that UTC, and the phase
advances at the J2 secular argument-of-latitude rate
du/dt = n [1 + 0.75 J2 (R_eq/a)² (8 cos² i − 2)] (draconic period).

```python
from sdrift_sky import SkyModel
sky = SkyModel()
p = sky.point(
    30, 10, "2026-10-17T00:10:00Z",
    model="total", orbit="noon-midnight", height=650,
    orbital_phase_deg=90, epoch="2026-10-17T00:00:00Z",
)
print(p["orbital_phase_deg"])  # about 126.847°
```

At 650 km the draconic period is 97.70 min (the Keplerian value is 97.58 min).
The reference UTC is a design assumption.

### 4.5 ELV and avoidance angles

The zenith angle θ is the angle between the satellite zenith and the line of
sight: 0° is the zenith, 90° horizontal, 180° the nadir. There is no minimum or
maximum zenith-angle filter.

For a spherical Earth of radius R = 6371 km:

```text
ELV = 180° − asin[R / (R + H)] − θ
```

ELV = 0° is tangent to the Earth's surface; positive values point above it.
At 650 km, `ELV ≈ 114.85° − θ`, so the default ELV ≥ 30° means θ ≤ 84.85°.

By default a direction must satisfy all three conditions:

- ELV ≥ 30°
- separation from the Sun ≥ 90°
- separation from the Moon ≥ 60°

Earth-occulted directions are also excluded. Turning off avoidance does not
remove Earth occultation, and a negative ELV limit does not return brightness
for directions blocked by the Earth.

The Moon direction is the one seen from the satellite (parallax up to about
1.1° in LEO). The Sun direction is geocentric (parallax from LEO < 0.003°).

## 5. Coordinate systems and direction input

### 5.1 `point(ra, dec, date, ...)`

**The input RA and DEC are always equatorial degrees.** Even with
`coords="galactic"`, the first two numbers of `point()` are not read as l, b.
RA is wrapped modulo 360°, and DEC must be within −90° to 90°.

### 5.2 `sky_map(..., coords=...)`

- `equatorial`: RA, DEC HEALPix map
- `galactic`: Galactic longitude l, latitude b HEALPix map

The map, Sun, Moon, zenith, and avoidance boundaries are shown in the same
coordinate system, and the FITS file stores the pixel grid of the selected
system. `point()` and `evaluate()` also accept `coords` as metadata for the
result, but it does not change the meaning of their input RA/DEC or direction
vectors.

### 5.3 Querying a single l, b direction

```python
from sdrift_sky import SkyModel, vector, coordinates
from sky_coordinates import to_equatorial

sky = SkyModel()
l, b = 120.0, 30.0
ra, dec = coordinates(to_equatorial(vector(l, b), "galactic"))
p = sky.point(ra, dec, "2026-10-17T00:00:00Z", model="total")
print(ra, dec, p["ab_mag_arcsec2"], p["visible"])
```

The conversion is a fixed rotation of the axes; it adds no parallax or
aberration correction.

## 6. Single-direction result: `point()`

The result is a Python `dict` of scalars.

| Key | Units, meaning |
|---|---|
| `ra_deg`, `dec_deg` | equatorial coordinates of the input direction, degrees |
| `utc`, `band`, `model`, `airglow_model`, `coords` | settings |
| `ab_mag_arcsec2` | AB mag/arcsec²; NaN for excluded directions |
| `intensity_mjy_sr` | total intensity of the selected components, MJy/sr |
| `zodiacal_mjy_sr`, `dgl_mjy_sr`, `airglow_mjy_sr` | components, MJy/sr |
| `photon_radiance` | airglow in photons/s/cm²/sr; 0 if airglow is not selected |
| `view_allowed` | ELV and Sun/Moon avoidance satisfied |
| `earth_occulted` | Earth occultation |
| `visible` | `view_allowed` and not Earth-occulted |
| `model_valid` | support range / cutoff check of the numerical model |
| `valid` | `visible and model_valid` |
| `display_ab_mag_arcsec2` | brightness if visible, otherwise `None` |
| `orbit`, `height_km`, `orbital_phase_deg` | orbit type, km, phase at this UTC in degrees |
| `sza_deg`, `zenith_deg`, `elv_deg`, `minimum_elv_deg` | satellite SZA, zenith angle, ELV, ELV limit, degrees |
| `sun_ra_dec_deg`, `moon_ra_dec_deg` | equatorial `(RA, DEC)` of the Sun (geocentric) and the Moon (seen from the satellite) |
| `sun_separation_deg`, `moon_separation_deg` | angle between the line of sight and the Sun/Moon, degrees |
| `sun_avoidance`, `moon_avoidance` | avoidance switches |
| `minimum_sun_separation_deg`, `minimum_moon_separation_deg` | avoidance angles, degrees |

`visible=True` can occur with `model_valid=False`: the direction is observable,
but the brightness is a partial result. `model_valid=True` also does not mean
that missing emission components are included. Do not interpret the model flag
for excluded directions.

## 7. All-sky map: `sky_map()`

```python
from pathlib import Path
from sdrift_sky import SkyModel

sky = SkyModel()
result = sky.sky_map(
    "2026-10-17T00:00:00Z",
    model="total", airglow_model="combined", band="VIS",
    nside=128, coords="galactic",
    orbit="noon-midnight", height=650, orbital_phase_deg=90,
    min_elv_deg=30, sun_separation=90, moon_separation=60,
)
Path("output/python").mkdir(parents=True, exist_ok=True)
sky.save_fits(result, "output/python/total_VIS_nside128_20261017_galactic.fits.gz")
fig = sky.plot(result, contour_step=0.2)
fig.savefig("output/python/total_galactic.png", dpi=160)
```

The settings are those of `point()`, plus `nside` and the map coordinate system.

### 7.1 Resolution and pixel order

`nside` must be a positive power of two; the default is 128. The GUI offers
16/32/64/128; the API and CLI accept any valid NSIDE.

```text
Number of pixels N = 12 × nside²
NSIDE 16  →   3,072
NSIDE 32  →  12,288
NSIDE 64  →  49,152
NSIDE 128 → 196,608
```

Arrays are in **RING HEALPix** order. A higher NSIDE costs more computation but
does not add information beyond the NSIDE 128 zodiacal/DGL input templates.
A Galactic map evaluates the actual line of sight of each Galactic pixel center.

### 7.2 Main keys of the map dictionary

Direction arrays have shape `(N,)`. The structure is the same for every `model`.

| Key | Content |
|---|---|
| `flux`, `intensity_mjy_sr` | the same total intensity array, MJy/sr |
| `flux_unit` | always `"MJy/sr"` |
| `mag`, `ab_mag_arcsec2` | the same AB mag/arcsec² array |
| `zodiacal`, `dgl`, `airglow` | component arrays, MJy/sr |
| `airglow_photon_radiance` | airglow photon radiance |
| `lit` | airglow photon radiance from the sunlit part of the path |
| `view`, `view_allowed` | the same ELV and Sun/Moon avoidance mask |
| `visible` | avoidance satisfied and not Earth-occulted |
| `earth`, `earth_occulted` | the same Earth-occultation mask |
| `model_valid` | support range check of the numerical model |
| `valid` | `visible & model_valid` |
| `theta`, `elv` | zenith angle and ELV, degrees |
| `hmin` | minimum altitude along the forward line of sight, km; 0 if Earth-occulted |
| `cutoff` | path crosses the airglow photoelectron cutoff region |
| `nside`, `coords`, `model`, `band`, `utc` | map settings |
| `airglow_model`, `orbit`, `height`, `epoch`, `min_elv_deg` | settings |
| `state` | satellite geometry; see below |
| `sun`, `moon` | unit vectors in equatorial axes, `(3,)` each; Moon seen from the satellite |
| `engine` | actual path: `numba`, `numpy`, or `map interpolation` |

Main keys of `state`:

| Key | Meaning |
|---|---|
| `zenith` | satellite zenith unit vector at this UTC, equatorial axes |
| `orbital_phase_deg` | phase after evolution from the epoch, degrees |
| `sza` | solar zenith angle at the satellite, degrees |
| `period` | draconic period, **seconds** |
| `inclination_deg`, `raan_deg`, `ltan_hours` | inclination, RAAN (mean equinox of date), LTAN |
| `limb_zenith_deg` | zenith angle of the direction tangent to the surface, degrees |
| `axes` | rotation matrix `(3,3)` from the integration reference frame to equatorial axes |

Pass the returned dictionary unchanged to `plot()` and `save_fits()`. Keys
starting with an underscore, such as `_engine_result`, are internal.

### 7.3 Statistics of observable directions

```python
import numpy as np
from sdrift_sky import SkyModel

sky = SkyModel()
r = sky.sky_map("2026-10-17T00:00:00Z", model="total", nside=16)
use = r["visible"] & r["model_valid"] & np.isfinite(r["mag"])
if use.any():
    print("Selected pixels:", np.count_nonzero(use))
    print("AB magnitude range:", r["mag"][use].min(), r["mag"][use].max())
else:
    print("No valid pixels satisfy the conditions.")
```

The mean of magnitudes differs from the magnitude of the mean intensity. For a
mean brightness, average `flux` and convert it with `ab_magnitude()`.

## 8. Direct evaluation of many lines of sight: `evaluate()`

`directions` are **unit vectors in equatorial axes**, shape `(3,)` or `(N,3)`.
Even a single vector input returns `(1,)` arrays. Do not pass RA, DEC numbers in
place of vectors.

```python
import numpy as np
from sdrift_sky import SkyModel, vector

sky = SkyModel()
ra = [0.0, 30.0, 120.0]
dec = [60.0, 10.0, -20.0]
vectors = np.stack([vector(a, d) for a, d in zip(ra, dec)])
r = sky.evaluate(
    vectors, "2026-10-17T00:00:00Z", model="total",
    orbit="noon-midnight", orbital_phase_deg=90,
)
print(r["mag"], r["visible"], r["model_valid"])
```

The result of `evaluate()` is not a HEALPix map, so do not pass it to `plot()`
or `save_fits()`. Use `sky_map()` for all-sky output.

## 9. Plotting and saving

### 9.1 `plot(result, target=None, contour_step=0.2, figure=None)`

Takes a `sky_map()` result and returns a Matplotlib `Figure`.

| Argument | Meaning |
|---|---|
| `result` | dictionary returned by `sky_map()` |
| `target` | optional equatorial `(RA, DEC)` to mark, degrees |
| `contour_step` | default 0.2 mag; 0.1 to 1.0 mag in steps of 0.1 |
| `figure` | existing Matplotlib Figure to draw into; a new one if omitted |

Save with `fig.savefig("map.png")` or `fig.savefig("map.pdf")`, or call
`matplotlib.pyplot.show()` to display it. In a headless environment, set
`MPLBACKEND=Agg` before starting Python, or call `matplotlib.use("Agg")` before
importing `matplotlib.pyplot`.

The color range is set from the valid brightness range; larger magnitudes are
darker. The Sun is a red star with a red dashed boundary, the Moon a yellow/orange
circle with an orange dashed boundary, the minimum ELV a green dashed line, and
the zenith a blue dot. Gray is Earth occultation, pale areas are excluded
directions, and hatching marks directions where the airglow result is partial.

### 9.2 `save_fits(result, filename)`

The file name must end with `.fits.gz`. Missing parent folders are created, and
an existing file of the same name is overwritten. API v2 writes the following
**14 columns** in the same order and units for every mode.

| Column | Units |
|---|---|
| `AB_MAG_ARCSEC2` | mag/arcsec² |
| `INTENSITY_MJY_SR` | MJy/sr |
| `ZODIACAL_MJY_SR` | MJy/sr |
| `DGL_MJY_SR` | MJy/sr |
| `AIRGLOW_MJY_SR` | MJy/sr |
| `AIRGLOW_PHOTON_RADIANCE` | photons/s/cm²/sr |
| `EARTH_OCCULTED` | 0/1 |
| `MODEL_VALID` | 0/1 |
| `MIN_ALTITUDE` | km |
| `ZENITH_ANGLE` | degree |
| `VIEW_ALLOWED` | 0/1 |
| `ELV` | degree |
| `VISIBLE` | 0/1 |
| `VALID` | 0/1 |

Columns are float32, so convert the 0/1 columns to bool after reading if needed.
`COORDSYS=C` is equatorial, `G` Galactic, and `ORDERING=RING`.

The header contains `APIVER=2`, `MODEL` (sky mode), `DATATYPE`, `AGMODEL`
(airglow model), `FLUXUNIT=MJy/sr`, `DATE-UTC` (with the trailing `Z`), band
limits, ELV and avoidance settings, and source hashes. Orbit keywords:

| Keyword | Meaning |
|---|---|
| `ORBIT` | `dawn-dusk` or `noon-midnight` |
| `ORBMODEL` | `SSO J2 CIRCULAR` |
| `SAT_H` | satellite height, km |
| `INCL` | inclination, deg |
| `RAAN` | right ascension of the ascending node, deg (mean equinox of date) |
| `LTAN` | local time of the ascending node, hours |
| `SAT_SZA` | solar zenith angle at the satellite, deg |
| `PHASE`, `PHASE0` | phase at this UTC; input phase |
| `PHREF`, `PHCONV` | `NORTH`, `NORTH_V2` |
| `PHMODE` | `MANUAL` or `EPOCH` |
| `EPOCH` | reference UTC; written only when `PHMODE='EPOCH'` |
| `PERIOD` | draconic period, s |
| `ELVMIN`, `ELVREF`, `EARTHR` | ELV limit, `SPHERICAL EARTH`, 6371 km |
| `MOONOBS` | `SATELLITE` (Moon seen from the satellite) |

Zodiacal output also records `REF-UTC` (frame epoch of the zodiacal template)
and `ROT-DEG`, the rotation of the template frame about the J2000 ecliptic
north pole.

```python
from astropy.io import fits

with fits.open("output/python/total_VIS_nside128_20261017_galactic.fits.gz") as hd:
    header = hd[1].header
    table = hd[1].data
    mag = table["AB_MAG_ARCSEC2"].ravel()
    visible = table["VISIBLE"].ravel().astype(bool)
    print(header["COORDSYS"], header["APIVER"], mag[visible].shape)
```

This is the **common Python API v2 schema**. Files saved by the web page use a
different column layout that depends on the model (see
[README_web.md](README_web.md)), and older files may differ as well, so select
columns by name rather than by number. The Sun, Moon, and zenith RA/DEC header
values are equatorial regardless of the map coordinate system.

### 9.3 NaN, infinity, and masks

- `NaN`: no value, e.g. directions excluded by avoidance/ELV or Earth occultation
- `+inf` magnitude: the computed intensity is 0
- `visible`: observing geometry check
- `model_valid`: support range check of the numerical model
- `valid`: both checks satisfied

To store NaN or infinity in JSON, replace them with `None`. The CLI `--json`
output for a single direction does this automatically.

## 10. Command line: the same options in two files

These two commands run **the same CLI**.

```sh
python3 sdrift_sky.py --help
python3 airglow_orbit.py --help
```

### 10.1 Single direction

```sh
python3 sdrift_sky.py --ra 0 --dec 60 \
  --date 2026-10-13T12:00:00Z --model zodiacal_dgl --band VIS
```

Give `--ra` and `--dec` together. A single-direction request alone does not
produce a map. With `--png`, `--fits`, etc. as well, the direction is queried
and then the map is computed. `--ra` and `--dec` are equatorial even with
`--coords galactic`.

### 10.2 Galactic total map

```sh
python3 sdrift_sky.py --model total --airglow-model combined \
  --date 2026-10-17T00:00:00Z --band VIS --coords galactic --nside 128 \
  --orbit noon-midnight --orbital-phase-deg 90 --height 650 --min-elv 30 \
  --contour-step 0.2 --fits output/python/total_galactic.fits.gz \
  --png output/python/total_galactic.png
```

### 10.3 Changing only the components

```sh
python3 sdrift_sky.py --model dgl --date 2026-10-17T00:00:00Z \
  --orbit noon-midnight --orbital-phase-deg 90 --png output/python/dgl.png
python3 sdrift_sky.py --model airglow --airglow-model glow \
  --date 2026-10-17T00:00:00Z --orbit noon-midnight --orbital-phase-deg 90 \
  --png output/python/airglow.png
```

### 10.4 Main CLI options and API equivalents

| CLI | Python API / meaning |
|---|---|
| `--model` | `model`; default `total` |
| `--airglow-model` | `airglow_model` |
| `--date`, `--utc` | `date`; CLI default is the current UTC |
| `--ra`, `--dec` | equatorial coordinates for `point()` |
| `--orbit` | `orbit` |
| `--height`, `--H` | `height` |
| `--orbital-phase-deg`, `--phase` | `orbital_phase_deg` |
| `--epoch` | `epoch`; fixed phase if omitted |
| `--min-elv`, `--elv-min` | `min_elv_deg` |
| `--sun-separation`, `--moon-separation` | avoidance angles |
| `--no-sun-avoidance`, `--no-moon-avoidance` | set the matching switch to False |
| `--band`, `--coords`, `--nside` | map settings of the same name |
| `--engine` | `engine` |
| `--data-dir` | `data_dir` of the constructor |
| `--contour-step` | `contour_step` of `plot()`; default 0.2 |
| `--show` | open a map window |
| `--png`, `--pdf`, `--fits` | output paths |
| `--json` | save the single-direction dict; needs `--ra`, `--dec` |
| `--gui` | open the common GUI |

If both the output options and the coordinates are omitted, the map is computed
and a status line is printed, but nothing is saved or shown. Add the outputs you
need. `step_km` is an API argument only; there is no CLI option.

## 11. GUI

```sh
python3 sdrift_sky.py --gui
python3 airglow_orbit.py --gui
```

Both commands open the same window. Initial settings can be passed:

```sh
python3 sdrift_sky.py --gui --model dgl --coords galactic \
  --orbit noon-midnight --orbital-phase-deg 90 --height 650 --min-elv 25
```

| Control | Meaning |
|---|---|
| UTC / Now / ±1 h / ±1 min | observing time, current time, time shift |
| Map | `model` |
| Coordinates | equatorial or Galactic display |
| Orbit / H (km) / Orbital phase from north φ (°) | orbit, height, phase |
| Band / Model | band and airglow emission set `airglow_model` |
| NSIDE | map resolution |
| Min ELV (deg) | minimum ELV |
| Evolve phase from epoch / Phase epoch UTC | advance the phase from a reference time |
| Sun/Moon avoidance / Min separation (deg) | avoidance switches and angles |
| Contour (mag) | contour step |
| Calculate map | compute the map with the current inputs |
| Save PDF / PNG / FITS.GZ | save the current result |

Changing Map or Coordinates, or pressing a time button, triggers a new
calculation. After editing other fields, press **Calculate map**. Changing the
contour step redraws the existing result. Moving the pointer over the map shows
coordinates in the selected system, brightness, zenith angle, and ELV; the
status line shows φ and the satellite SZA. Calculations run in a worker thread;
the first Numba compilation can take longer. `--ra` and `--dec` are CLI query
options; the GUI has no target input field.

## 12. Faster and repeated calculation

`engine="auto"` uses parallel Numba airglow integration when Numba is available
and the NumPy path otherwise. `engine="numpy"` selects NumPy explicitly.
`engine="numba"` without Numba is an error for the `airglow` and `total` modes;
the zodiacal/DGL modes interpolate the input maps and do no line-of-sight
integration, so the engine setting does not matter there.

Reusing one `SkyModel` avoids repeated data reads and hash checks. Checking
inputs with a small NSIDE before the final NSIDE 128 run is convenient. For many
directions, one `evaluate()` call is faster than repeated `point()` calls.

## 13. Coordinate helper functions

| Function | Interface |
|---|---|
| `sdrift_sky.vector(ra, dec)` | two angles in degrees to a unit vector `(3,)` |
| `sdrift_sky.coordinates(v)` | one unit vector to `(longitude, latitude)` in degrees |
| `sky_coordinates.from_equatorial(v, coords)` | rotate equatorial vectors `(…,3)` to the selected system |
| `sky_coordinates.to_equatorial(v, coords)` | rotate vectors of the selected system to equatorial axes |
| `sdrift_sky.utc_time(date)` | returns `(Unix seconds, UTC string)` |
| `sdrift_sky.ab_magnitude(flux)` | MJy/sr to AB mag/arcsec² |

`vector()` is the generic longitude/latitude formula, so it can also convert
l, b as in Section 5.3. `coordinates()` does not know the axes of the vector, so
a vector in Galactic axes returns l, b.

## 14. Migrating from earlier code

| Earlier | Common API v2 |
|---|---|
| `SkyModel` and `AirglowOrbitModel` used differently | use `from sdrift_sky import SkyModel` |
| `kind="total"` | `model="total"` |
| `kind="airglow", model="glow"` | `model="airglow", airglow_model="glow"` |
| `--kind total --model combined` | `--model total --airglow-model combined` |
| only airglow `flux` in photon units | every `flux` in MJy/sr; photons in separate fields |
| earlier SkyModel returned brightness in excluded directions | excluded directions are NaN in every mode |
| earlier SkyModel had no orbit conditions | default orbit and ELV conditions in every mode |
| FITS columns differed by mode | common Python API v2 uses the same 14 columns |
| ideal solar-following orbit planes, `PHCONV=NORTH_V1` | J2 sun-synchronous orbits, `PHCONV=NORTH_V2` |
| Keplerian period | draconic period |
| geocentric Moon | Moon seen from the satellite |

`from airglow_orbit import AirglowOrbitModel` is a compatibility name for
`SkyModel` and uses **the same API and units**. Both script CLIs are identical.
The compatibility paths for `kind=`, `phase=`, `--kind`, `--phase`, and the old
GLOW model names remain, but use the recommended names in new code. Giving
conflicting `model` and `kind` is an error.

The earlier example `sky.point(0, 60, date, band="VIS", model="zodiacal_dgl")`
still runs, but follows the new orbit/ELV defaults and NaN policy. When
comparing with earlier numbers, check that defaults, avoidance conditions,
phase convention, orbit model, and units match. The orbit plane now follows the
J2 sun-synchronous geometry, so maps made with NORTH_V1 are not reproduced
exactly at the same φ. If an older noon-midnight phase u was measured from noon,
the north-referenced φ is `(u − 90°) mod 360°`.

## 15. Troubleshooting

| Symptom | What to check |
|---|---|
| `No module named ...` | required Python files and installed requirements |
| `manifest.json` not found | complete `docs/data/` or the `data_dir` path |
| SHA-256 mismatch | code/data versions and file damage; restore the data files together |
| `Lunar data must be geocentric vectors in km` | `docs/data/moon_gcrs.f32.gz` is from an older export; update `docs/data/` |
| `Unknown sky model` etc. | roles and spelling of `model` and `airglow_model` |
| Coordinate errors | RA/DEC in degrees, DEC within −90 to 90, both values given |
| NaN brightness | `visible`, ELV, Sun/Moon separation, Earth occultation |
| Infinite brightness | computed intensity is 0 |
| Entire map empty | overlap of the observable region at this phase with the ELV and avoidance limits |
| Hatching | airglow support range or photoelectron cutoff; not an observing restriction |
| GUI does not open | graphical desktop, Tkinter, Matplotlib backend |
| Numba errors | Python/NumPy/Numba compatibility; use `engine="numpy"` or `--engine numpy` |
| Column numbers differ from older files | check `APIVER`, column names and units |

## 16. Scope of the physical models

The date sets the sky positions, the Sun and Moon directions, the orbit plane,
and (in epoch mode) the orbital phase. Airglow emission rates come from a
precomputed table for mean solar activity; this is not a forecast model with
real-time F10.7/Ap or all atmospheric variations.

The zodiacal light is the SKYSURF model (O'Brien et al. 2025, a Kelsall-type
model). Each band template stores annual Fourier harmonics (order 3) in the
Sun's J2000 ecliptic longitude L, in a frame co-rotating with the Sun; for a
given date the harmonics are combined at L and the map is sampled in the
rotated frame. This includes the annual variation due to Earth's position
relative to the tilted, offset dust cloud and to the Earth-Sun distance. Below
0.25 μm the SKYSURF albedo and phase-function fits are extrapolations of the
model calibration (this affects the 0.20–0.25 μm part of UV1).

DGL uses the C2022_all correlation spectrum and the self-absorption correction
assumptions; the extension to short UV wavelengths is an assumption, not an
observation. Airglow still has missing components, no UV absorption, no
dawn/dusk chemistry history, and applies the PALACE night layers globally.
The Moon direction includes the satellite parallax; scattered moonlight is not
included. The orbit is circular with J2 secular rates only, and Earth's shadow
is a cylinder of radius R without penumbra, atmospheric extension, or
oblateness.

Further reading:

- [Building the datasets from the sources](README_first.md)
- [Orbit, ELV geometry and background](README_airglow_orbit.md)
- [Web usage](README_web.md)
- [Dawn-dusk phase figure](docs/assets/geometry_phase_dd.svg)
- [Noon-midnight phase figure](docs/assets/geometry_phase_nm.svg)
- [Satellite, Earth limb, and ELV geometry figure](docs/assets/geometry_limb_elv.svg)
- [GLOW license](docs/data/glowlicense.txt)

# S-DRIFT Sky Background — web and Python

To rebuild the data from the original sources, see the
[dataset generation guide](README_first.md).

`docs/index.html` computes zodiacal light, DGL, and airglow in the browser.
Users do not need Python, healpy, or NumPy. Current Chrome / Edge / Firefox /
Safari are supported; the page uses Web Workers, CompressionStream,
DecompressionStream, and Web Crypto. It works without an external CDN or a
calculation server; all data are read from `docs/data/`. One Web Worker stays alive
for the session and keeps the decompressed, SHA-256-checked arrays, so changing
the date or phase does not reload them; a newer request or Cancel stops the
running calculation.

## GitHub Pages

In **Settings → Pages → Build and deployment** of this repository, choose
**Deploy from a branch → main (or the actual default branch) → /docs → Save**.
Address after deployment: https://seoncafe.github.io/SDRIFT_sky/

Only `docs/` and this document are needed for the web page. `sdrift_sky.py` and
`requirements-python.txt` are also provided for Python users. The original
FITS/NPZ files, the Python GUI, the GLOW/PALACE sources, PDF results, tests, and
generation scripts are not needed for the web page. Include `docs/.nojekyll` as
well. All asset paths are relative.

Local preview:

```sh
python3 -m http.server 8765 --directory docs
```

Open http://localhost:8765/. Opening the HTML directly with `file://` may block
reading the data and modules; use HTTPS or localhost.

## Python: computing RA, DEC, UTC and displaying maps

For function arguments, return values, units, CLI, and GUI details, see the
[Python interface guide](README_python_interface.md). `SkyModel` is the common
entry point for all modes.

Python 3.10 or later is required. Download the repository and run it; the
repository's `docs/data/` is used without the original FITS files or a
GLOW/PALACE installation. No internet connection is needed at run time.

```sh
python3 -m pip install -r requirements-python.txt
python3 sdrift_sky.py --ra 0 --dec 60 --date 2026-10-13T12:00:00Z --band VIS --model zodiacal_dgl --show
```

RA and DEC are in **degrees**, and the date is in **UTC**. `--model` is
`zodiacal`, `dgl`, `zodiacal_dgl`, `airglow`, or `total`; the default is the
three-component sum `total`. The airglow emission set is chosen separately with
`--airglow-model combined/glow`. `--band` is UV1/UV2/UV3/VIS, default VIS.
Dates must be within 1900–2100. Omit `--show` to compute a single direction
only. The target is marked with a blue X on the map.

Saving files without a window:

```sh
python3 sdrift_sky.py --ra 0 --dec 60 --date 2026-10-13T12:00:00Z \
  --band VIS --nside 128 --contour-step 0.2 \
  --pdf sky.pdf --png sky.png --fits sky.fits.gz --json point.json
```

All modes share the orbit, height, phase, minimum ELV, the Sun/Moon positions,
the red/orange dashed separation boundaries, and the avoidance masks. The
default Sun separation is 90° and the Moon separation 60°; change them with
`--sun-separation 100 --moon-separation 70`, or turn them off with
`--no-sun-avoidance` and `--no-moon-avoidance`. The contour step defaults to
0.2 mag; `--contour-step` selects 0.1–1.0 mag. FITS.gz stores excluded
directions as NaN. Python FITS files use the same 14 columns in every mode,
including each component (see the Python interface guide).

Direction queries and maps apply the same observing conditions. Excluded
directions have NaN brightness (null in JSON); check `visible`, `valid`, and the
actual Sun/Moon separations. `flux`/`intensity_mjy_sr` are in MJy/sr in every
mode; the airglow photon radiance is kept in `photon_radiance` for direction
queries and in `airglow_photon_radiance` for maps. Maps are computed at HEALPix
pixel centers, so a value read with the web pointer and a direct query at the
same coordinates can refer to directions up to a pixel apart. Brightness is
shown in AB mag / arcsec².

Calling from Python:

```python
from sdrift_sky import SkyModel
sky = SkyModel()
point = sky.point(0, 60, "2026-10-13T12:00:00Z", band="VIS", model="zodiacal_dgl")
print(point)
result = sky.sky_map("2026-10-13T12:00:00Z", nside=128, band="VIS", model="zodiacal_dgl")
figure = sky.plot(result, target=(0, 60), contour_step=0.2)
figure.savefig("sky.pdf", bbox_inches="tight")
```

`--show` needs a GUI backend and a graphical desktop; on a server use
`--pdf`/`--png`. The Sun direction is the stored geocentric GCRS data shared
with the web page. The Moon direction is computed from the stored geocentric
Moon vectors (km) as seen from the satellite, identically in Python and in the
browser; scattered moonlight is not included. DGL follows the α ≈ α′
normalization assumption described below, and the UV extension of the
correlation is not an observation.

## Using the web page

- Models: Zodiacal light / DGL / Zodiacal light + DGL / Airglow / Zodiacal light
  + DGL + Airglow. Bands: UV1 (200–400 nm), UV2 (250–400 nm), UV3 (300–400 nm),
  VIS (400–700 nm).
- The default model is Zodiacal light + DGL + Airglow, starting at the current
  UTC. The ▲/▼ buttons next to each date/time field and the keyboard ↑/↓ keys
  change the value and carry over field boundaries (month to year, day to
  month, second to minute, etc.). Year and month changes keep the day within
  the target month (January 31 plus one month is February 28 or 29).
- Every model uses the same orbit inputs: orbit (dawn-dusk / noon-midnight),
  satellite height H (default 650 km), phase φ from the northernmost point
  (default 0°), minimum ELV (default 20°), and an optional phase evolution from
  a reference UTC. The panel shows the sun-synchronous inclination i for the
  entered height (97.96° at 650 km, 96.65° at 300 km, 99.45° at 1000 km); it
  updates as H changes. See [README_airglow_orbit.md](README_airglow_orbit.md)
  for the orbit model and the phase convention.
- With Sun avoidance on, the Sun separation boundary is a red dashed line. The
  Sun is always marked with a red star and the label Sun, even with Sun
  avoidance off. The default Sun separation is 90°, in 5° steps.
- The Moon position and its avoidance boundary are shown in every model. Moon
  avoidance is on by default, with a minimum separation of 60° in 5° steps. With
  both Sun and Moon avoidance on, only directions satisfying both are shown and
  saved.
- The Moon is computed with the Astropy builtin ephemeris as geocentric GCRS
  vectors in km every 3 hours, interpolated, and then viewed from the satellite
  (parallax up to about 1.1° in LEO). 1900–2100 are supported; scattered
  moonlight is not included. FITS headers record MOONMASK, MOONMIN, MOONRA,
  MOONDEC, and MOONOBS=SATELLITE.
- DGL and Zodiacal light + DGL use the same date/time and Sun avoidance settings
  as the zodiacal light. The DGL distribution is fixed in RA, DEC; the two
  intensities are added and then converted to AB magnitude. The summed FITS also
  contains each component intensity.
- Airglow depends on the date through the orbit geometry: the Sun direction,
  the orbit plane, and the satellite position (orbit, H, and φ) set the local
  SZA along every line of sight. The emission set is GLOW + PALACE night layers
  by default. The emission rates are mean-solar-activity tables and do not change
  with the date.
- NSIDE defaults to 128; choose 32 or 64 on slow devices. Changing an input
  starts a new calculation automatically, and the stop button cancels it.
- RA increases to the left. The contour step defaults to 0.2 mag and can be
  set from 0.1 to 1.0 mag in 0.1 mag steps. A change is applied to the current
  map immediately and to saved PDF/PNG files.
- Each contour level is labeled. Line widths increase from 0.6 to 2.2 px from
  the lowest to the highest level, identically on screen and in PDF/PNG output.
  Larger AB magnitudes are darker. The pointer shows coordinates, the value,
  the zenith angle θ, and ELV.
- The coordinate system can be Equatorial (RA, DEC) or Galactic (l, b); FITS
  files are saved in the selected system (see below).
- PDF / print: choose PDF as the destination in the print dialog. The map title,
  legend, and model license notes are included. PNG can also be saved.
- File names follow the pattern
  `{model}_{band}_nside{N}_{YYYYMMDD_HHMMSS}_{orbit}_H{H}_phi{φ:.3f}_elv{ELV}[_{airglow set}]_{coords}[_moon{angle}][_anti-sun]`,
  for example
  `zodiacal_VIS_nside128_20261008_123000_dawn-dusk_H650_phi0.000_elv20_equatorial_moon60_anti-sun.fits.gz`.
  The time is UTC; the airglow emission set (`combined` or `glow`) appears only
  for the Airglow and total models; `_moon…` and `_anti-sun` appear when Moon
  and Sun avoidance are on.

## Web FITS layout

The browser writes a RING HEALPix BINTABLE with `COORDSYS=C` (equatorial) or
`G` (Galactic), `WEBVER=2`, and `CREATOR='S-DRIFT HTML'`. Its columns depend on
the model and differ from the fixed 14-column layout of Python API v2:

| Model | Columns |
|---|---|
| `zodiacal`, `dgl` | AB_MAG_ARCSEC2, INTENSITY_MJY_SR, EARTH_OCCULTED, MODEL_VALID, MIN_ALTITUDE, ZENITH_ANGLE, VIEW_ALLOWED, ELV |
| `zodiacal_dgl` | AB_MAG_ARCSEC2, INTENSITY_MJY_SR, ZODIACAL_MJY_SR, DGL_MJY_SR, EARTH_OCCULTED, MODEL_VALID, MIN_ALTITUDE, ZENITH_ANGLE, VIEW_ALLOWED, ELV |
| `airglow` | AB_MAG_ARCSEC2, PHOTON_RADIANCE, SUNLIT_RADIANCE, EARTH_OCCULTED, MODEL_VALID, MIN_ALTITUDE, PHOTO_CUTOFF, ZENITH_ANGLE, VIEW_ALLOWED, ELV |
| `total` | AB_MAG_ARCSEC2, INTENSITY_MJY_SR, SUNLIT_RADIANCE, EARTH_OCCULTED, MODEL_VALID, MIN_ALTITUDE, PHOTO_CUTOFF, ZENITH_ANGLE, VIEW_ALLOWED, ELV, ZODIACAL_MJY_SR, DGL_MJY_SR, AIRGLOW_MJY_SR, AIRGLOW_PHOTON_RADIANCE |

PHOTON_RADIANCE, SUNLIT_RADIANCE, and AIRGLOW_PHOTON_RADIANCE are in
photons/s/cm²/sr; SUNLIT_RADIANCE is the airglow part from the sunlit section of
the path. EARTH_OCCULTED and MIN_ALTITUDE are filled for every pixel, as in
Python. The orbit keywords are the same as in Python: ORBIT,
ORBMODEL='SSO J2 CIRCULAR', SAT_H, INCL, RAAN (mean equinox of date), LTAN,
SAT_SZA, PHASE, PHASE0, PHREF=NORTH, PHCONV=NORTH_V2, PHMODE, PERIOD (draconic
period, s), ELVMIN, and ZENRA/ZENDEC. EPOCH is written only when
PHMODE='EPOCH'. DATE-UTC keeps the trailing `Z`. Read columns by name.

## What is computed

AB mag/arcsec² for a constant **photon-counting** throughput.

The zodiacal light is the SKYSURF model (O'Brien et al. 2025, a Kelsall-type
model). Each band template stores annual Fourier harmonics in the Sun's J2000
ecliptic longitude L, in a frame co-rotating with the Sun:
`I = A0 + Σ_{k=1..3} [Ck cos kL + Sk sin kL]`. The browser (`zodiacalTemplate`
in `sky-core.mjs`) and Python (`zodiacal_template` in `sdrift_sky.py`) combine
the harmonics for the date and then sample the map rotated about the J2000
ecliptic pole by L − L_ref; the template NSIDE is read from the manifest. This
includes the annual variation from Earth's position relative to the tilted,
offset dust cloud and the Earth-Sun distance. Below 0.25 μm the SKYSURF albedo
and phase-function fits are extrapolations of the model calibration. The Sun
vectors are the same Astropy get_sun vectors in J2000 axes as in Python, stored
every 12 hours and interpolated. Supported dates are UTC 1900–2100. Sun
avoidance keeps directions with Sun separation ≥ the given value; the default is
90°.

DGL uses IRIS 100 μm, SFD E(B−V), and Fitzpatrick 1999 with R_V = 3.1. The Å
values of the C2022_all table are converted to μm and interpolated, and the
original correction factor 2.1 is kept. Below 372.5 nm the correlation is
extended as a constant 0.127596. The self-absorption formula of the paper's
Eqs. (3)–(6), `I_nu = alpha * beta * I100 * lambda_um/100`, is applied.
However, C2022_all is the alpha_prime table of Figure 3, and the original alpha
and the effective sample beta are not available, so **beta_ref = 1 (optically
thin sample) and alpha ≈ alpha_prime** are assumed. The self-absorption law
follows the paper, but the normalization is an approximation. The correction
factor 2.1 is kept, and the earlier inverse-beta maps are not used. The UV
extension and the extrapolation to the Galactic plane / high optical depth are
assumptions. The wavelength integral assumes constant photon counting with
steps of at most 0.5 nm, and the DGL does not rotate with the date. The color
range is set from the minimum and maximum valid magnitudes after masking,
rounded outward to 0.5 mag. FITS records BETAOP=MULTIPLY, ALPHAPP=T, and
BETAREF=1. Source: https://arxiv.org/abs/2201.01378 .

Airglow uses the emission rates q(h, SZA, component) summed into two channels for
each band: day/all and night-only. Without absorption and with a constant
throughput, this sum is equivalent to integrating the selected components
separately. The browser interpolates bilinearly in altitude and local SZA and
integrates each HEALPix line of sight with 2 km midpoint steps. It does not
re-interpolate the Python H/SZA integral LUT, so it has none of that LUT's
additional error; float32 storage rounding remains.

Each line of sight is integrated in a reference frame with the Sun on the +X
axis and the satellite at `(R+H) × (cos SZA, 0, sin SZA)`, where SZA is the
satellite SZA of the selected orbit, date, and phase. The frame is rotated to
the actual sky (`orbitState(...).axes` in `sky-core.mjs`), so the map is in real
RA/DEC (or l, b) for that date and orbit. R = 6371 km, and the table top is
1500 km. The PALACE night components are used only in Earth's shadow (a
cylinder of radius R, without penumbra, atmospheric extension, or oblateness).
Time-dependent chemistry and twilight history are not included.

The airglow band components are **incomplete**. In particular, the 200–300 nm
molecular and continuum spectra of UV1/UV2, UV line-of-sight absorption, the H
geocorona, and some high vibrational levels / rotational spectra are missing.
PALACE applies the Paranal nighttime column inside Earth's whole shadow and
assumes an altitude distribution. The GLOW background is F10.7 = 100, Ap = 4;
for each SZA χ it is the mean over four seasons and the two local solar times
12 h ∓ χ/15 h (equinox, equator relation). It is not a global observed
climatology. Gray marks Earth occultation; hatching marks directions where the
model is incomplete because of the minimum supported altitude or the GLOW
photoelectron cutoff (SZA ≥ 1.85 rad). Zero radiance is +∞ mag and is excluded
from the finite color range and contours.

Sources, components, inputs, and SHA-256 hashes are in `docs/data/manifest.json`.
PALACE: https://gmd.copernicus.org/articles/18/4353/2025/ (CC BY 4.0).
The GLOW and zodiacal licenses are included in `docs/data/`.

This software is based on calculations performed using the GLOW model.
Use of the GLOW model is governed by the Open Source Academic Research License
Agreement contained in [glowlicense.txt](docs/data/glowlicense.txt).
Browser routines modified 2026-10-20.

## Airglow maps for the two sun-synchronous orbits

In every model, choose dawn-dusk / noon-midnight, UTC, height, phase, and
minimum ELV. Both orbits show the Sun and Moon positions and boundaries; the
default Sun/Moon separations are 90°/60°. Start the Python GUI with
`python3 airglow_orbit.py --gui`. FITS.gz, PDF, PNG output and CLI calculation
are supported. For the orbit model (circular sun-synchronous orbit with J2,
LTAN 18:00 / 12:00), phase evolution, examples, and limitations, see
[README_airglow_orbit.md](README_airglow_orbit.md).

The Zodiacal light + DGL + Airglow mode adds the intensities of the three
components and uses the same orbit, angles, ELV, and Sun/Moon avoidance as
Airglow. Select it in the web model list or run
`python3 airglow_orbit.py --gui --model total`.
[Usage](README_airglow_orbit.md#zodiacal-light--dgl--airglow-sum).

The default mode of both the web page and the Python orbit viewer is
**Zodiacal light + DGL + Airglow**.

The coordinate selector switches between **Equatorial (RA, DEC)** and
**Galactic (l, b)** maps. The Python orbit GUI starts in Galactic coordinates
with `--coords galactic`, and both Python CLIs and
`sky_map(..., coords='galactic')` support it. FITS COORDSYS and the actual pixel
grid follow the selected system. For Python, keep the common transformation
module `sky_coordinates.py` in the same folder.

Zodiacal, DGL, and Zodiacal + DGL also accept all orbit, phase, height, and ELV
settings. The Python GUI/CLI supports `airglow_orbit.py --model zodiacal`,
`--model dgl`, and `--model zodiacal_dgl`, applying the common observing
conditions to the brightness of the selected components.

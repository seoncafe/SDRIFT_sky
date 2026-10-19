# Date- and orbit-dependent airglow maps

For the early design phase without a TLE, two circular sun-synchronous orbits
are compared. Python and the web page use the same mean emission tables and
Sun/Moon position data in `docs/data/`; GLOW or PALACE do not need to be
installed. Output coordinates are RA, DEC in geocentric GCRS axes, as for the
zodiacal/DGL maps. The default height is 650 km and the default NSIDE is 128.

For the common Python API, CLI, and GUI, see the
[Python interface guide](README_python_interface.md). `SkyModel` and the
compatibility name `AirglowOrbitModel` are the same class.

## Running

```sh
python3 -m pip install -r requirements-python.txt
# Optional: a Numba compatible with the installed Python/NumPy enables parallel integration
python3 -m pip install numba
python3 airglow_orbit.py --gui
```

`python3 airglow_orbit_gui.py` opens the same window. In the development folder,
which also contains the source calculation code, `python3 airglow_gui.py` opens
this viewer too, and `python3 airglow_gui.py --reference` starts the earlier
fixed-H/SZA reference viewer. The GitHub release contains only the new viewer.
The GUI starts at the current UTC in the Zodiacal + DGL + Airglow mode. Choose
the orbit, UTC, H, phase φ from the northernmost point, band, and minimum ELV,
then press Calculate map. Contours range from 0.1 to 1.0 mag; the default is
0.2 mag. Tkinter and a graphical desktop are required. Check with
`python3 -m tkinter`; if Tkinter is missing, install the Tcl/Tk package matching
your Python distribution. With Homebrew Python on macOS, use the `python-tk@3.x`
matching that Python version. Do not install `tkinter` with pip.

On the web, select the **Airglow** model in `docs/index.html`.

```sh
python3 -m http.server 8765 --directory docs
```

Open http://localhost:8765/. GitHub Pages serves the same interface. The lines
of sight are integrated in a Web Worker, and the calculation can be stopped.
Python integrates in parallel when Numba is available and uses NumPy otherwise;
`--engine numpy` selects NumPy explicitly.

## Zenith angle, Sun and Moon avoidance

The zenith angle θ is measured from the satellite zenith (the direction away
from Earth's center): 0° is the zenith, 90° horizontal, 180° toward Earth's
center. The zenith angle is shown in the pointer readout and stored in the FITS
file for reference; there are no minimum/maximum zenith-angle inputs or filters.
The observable region of every line of sight is set by the minimum ELV, Sun and
Moon avoidance, and Earth occultation. Downward lines of sight not blocked by
Earth also integrate the atmosphere below the satellite.

**The default observing conditions are Sun separation ≥ 90°, Moon separation
≥ 60°, and ELV ≥ 20°.** ELV (Earth limb elevation) is the angle of the line of
sight above the limb of the spherical Earth. For height H and zenith angle θ,
`ELV = 180° − asin(R/(R+H)) − θ` with R = 6371 km. At 650 km,
`ELV ≈ 114.85° − θ`. Set the **minimum ELV** in the GUI/web, with
`--min-elv 10` in the CLI, or `min_elv_deg=10` in the Python API. The default
is 20.0°, and the allowed range is −90° to 180°. It is applied together with
Sun and Moon avoidance. Even with a negative value, directions blocked by Earth
are not computed. ELV ≥ 0° does not mean that the line of sight avoids every
atmospheric layer of the model. The minimum ELV boundary is a green dashed line,
and the pointer readout shows ELV and the zenith angle.

Python and the web show the Sun and the Moon; the Sun boundary is a red dashed
line and the Moon boundary an orange dashed line. The avoidance angles can be
changed or turned off. The zenith is shown in blue. Directions outside the
avoidance/ELV conditions are pale gray, Earth occultation is gray, and outside
the model range or beyond the photoelectron cutoff is hatched. Hatched
brightness is a partial result.

The Moon direction is the one seen from the satellite (parallax up to about
1.1° in LEO), so Moon avoidance and the Moon marker use that direction. The Sun
direction is geocentric (parallax from LEO < 0.003°).

## Orbit model

Both orbits are circular, with radius R + H (R = 6371 km, spherical Earth), and
include only the J2 secular rates.

- **Inclination.** From the J2 sun-synchronous condition
  `dΩ/dt = −1.5 n J2 (R_eq/a)² cos i = 360° per tropical year`
  (J2 = 1.08262668e-3, R_eq = 6378.137 km, μ = 398600.4418 km³/s²).
  At H = 650 km, i = 97.96°.
- **Node.** The ascending node is at a fixed local time in **mean** solar time:
  LTAN 18:00 for dawn-dusk and 12:00 for noon-midnight, so
  `RAAN = mean-Sun RA (of date) + 15° × (LTAN − 12 h)`. The plane follows the
  mean Sun; the true Sun differs by the equation of time and by its
  declination.
- **Pole.** Earth's pole is the mean pole of date (IAU 1976 precession; agrees
  with ERFA to < 0.3 arcsec over 1900–2100). Nutation (< 0.005°) and the frame
  bias are neglected.

Because the plane is inclined and follows the mean Sun, **the satellite SZA
varies around the orbit**. Dawn-dusk at 650 km: about 82–98° at the equinoxes,
74.5–105.5° at the December solstice, and 58.6–121.4° at the June solstice.
Near the June solstice the satellite therefore enters Earth's shadow (SZA
> 114.85° at 650 km). Noon-midnight SZA spans about 2–178°. The current SZA is
shown in the map title and stored as `SAT_SZA`.

Earth's shadow (used for the night layers and for the sunlit part of each path)
is a cylinder of radius R behind Earth. This is a deliberate approximation: no
penumbra (about 60 km wide at LEO distances), no atmospheric extension of the
shadow (tens of km), and no oblateness.

## Orbital phase φ (convention NORTH_V2)

The position on the orbit is not determined by the date alone, so it is set by
the phase φ, measured from the northernmost point of the orbit and advancing in
the direction of motion.

| φ | Dawn-dusk | Noon-midnight |
|---|---|---|
| 0° | northernmost point | northernmost point |
| 90° | descending node, dawn side | descending node, midnight side |
| 180° | southernmost point | southernmost point |
| 270° | ascending node, dusk side | ascending node, noon side |

The northernmost point is the maximum-declination point of the orbit (about 82°
for i = 97.96°), not a position above the geographic pole. If an older
noon-midnight phase u was measured from noon, φ = (u − 90°) mod 360°.
**The default phase is 0° (northernmost) for both orbits.** For the center of
the noon-midnight night side choose φ = 90°. At noon (φ = 270°) the observable
region is limited by Sun avoidance and the ELV limit.

By default the selected phase is kept when the date changes; the orbit plane,
the Sun, and the Moon still move with the date. To move the satellite with time,
select **evolve phase from the reference UTC** in the web page / **Evolve phase
from epoch** in the GUI, or add `--epoch` in the CLI. The phase at the reference
UTC is the given `--orbital-phase-deg`, and it advances at the J2 secular
argument-of-latitude rate

```text
du/dt = n [1 + 0.75 J2 (R_eq/a)^2 (8 cos^2 i − 2)],   n = sqrt(mu / a^3)
```

(draconic period). At 650 km the period is 97.70 min (the Keplerian value is
97.58 min). The Python API variable is `orbital_phase_deg`, in degrees. The
older `--phase` and API `phase=` are accepted as aliases and are interpreted in
the NORTH_V2 convention. Existing files are not modified. The default reference
UTC 2026-10-17 00:00:00 in the examples is an arbitrary design reference. In
epoch mode a date change moves the Sun, the Moon, the orbit plane, and the
orbital phase.

This is a design comparison model, not orbit propagation: there is no TLE/SGP4,
drag, or higher-order gravity term. References for orbits and conventions:
[NASA orbit overview](https://science.nasa.gov/earth/earth-observatory/catalog-of-earth-satellite-orbits/),
[Astropy satellite coordinates / TLE handling](https://docs.astropy.org/en/stable/coordinates/satellites.html).

## CLI and Python API

```sh
# Dawn-dusk, northernmost point at the selected UTC, default ELV and Sun/Moon avoidance
python3 airglow_orbit.py --model airglow --date 2026-10-17T00:00:00Z \
  --orbit dawn-dusk --height 650 --orbital-phase-deg 0 --band VIS --nside 128 \
  --fits airglow_VIS_nside128_20261017_dawn-dusk.fits.gz \
  --png airglow_VIS_nside128_20261017_dawn-dusk.png

# Noon-midnight, midnight side; save a PDF or use --show for a window
python3 airglow_orbit.py --model airglow --date 2026-10-17T00:00:00Z \
  --orbit noon-midnight --orbital-phase-deg 90 --height 650 --band UV3 \
  --pdf airglow_UV3_nside128_20261017_noon-midnight.pdf

# Phase evolved from the epoch, single RA/DEC query
python3 airglow_orbit.py --model airglow --date 2026-10-17T00:10:00Z \
  --epoch 2026-10-17T00:00:00Z --orbital-phase-deg 90 --orbit noon-midnight \
  --ra 30 --dec 10 --band VIS
```

`--sun-separation 100 --moon-separation 70` changes the avoidance angles;
`--no-sun-avoidance` and `--no-moon-avoidance` turn them off. Single-direction
queries apply the same conditions. Excluded directions return NaN brightness,
together with `view_allowed`, `earth_occulted`, `model_valid`, the zenith
angle, and the separations.

```python
from airglow_orbit import AirglowOrbitModel
sky = AirglowOrbitModel()
r = sky.sky_map('2026-10-17T00:00:00Z', model='airglow', orbit='noon-midnight',
                orbital_phase_deg=90, height=650, band='VIS', nside=128, min_elv_deg=20)
fig = sky.plot(r, contour_step=0.2)
fig.savefig('map.pdf')
sky.save_fits(r, 'map.fits.gz')
```

The Python FITS.gz file is a RING HEALPix table in the selected coordinate
system (`COORDSYS=C` for equatorial, `G` with `--coords galactic`). It uses the
fixed 14-column API v2 layout in every mode: AB_MAG_ARCSEC2, INTENSITY_MJY_SR,
ZODIACAL_MJY_SR, DGL_MJY_SR, AIRGLOW_MJY_SR, AIRGLOW_PHOTON_RADIANCE,
EARTH_OCCULTED, MODEL_VALID, MIN_ALTITUDE, ZENITH_ANGLE, VIEW_ALLOWED,
ELV (deg), VISIBLE, VALID (details in the
[Python interface guide](README_python_interface.md#92-save_fitsresult-filename)).
Brightness is NaN in excluded regions. VIEW_ALLOWED checks the minimum ELV and
Sun/Moon avoidance; Earth occultation and model validity are separate columns.
Excluded directions are not integrated, so their other diagnostic columns are
not meaningful.

The header records the UTC (`DATE-UTC`), orbit (`ORBIT`,
`ORBMODEL='SSO J2 CIRCULAR'`, `INCL`, `RAAN` in the mean equinox of date,
`LTAN` in hours), phase (`PHASE`, `PHASE0`, `PHREF=NORTH`, `PHCONV=NORTH_V2`,
`PHMODE`), the draconic `PERIOD` in seconds, height, satellite SZA, avoidance
angles, minimum ELV (`ELVMIN`), the RA/DEC of the Sun, Moon, and zenith, the
SHA-256 of the emission table, and `MOONOBS=SATELLITE`. `EPOCH` is written only
in epoch mode (`PHMODE='EPOCH'`). File names written by the GUI and the web page
mark the orbital phase with `_phi...`.

The web page saves a different, model-dependent column layout (for airglow it
includes PHOTON_RADIANCE, SUNLIT_RADIANCE, and PHOTO_CUTOFF); see
[README_web.md](README_web.md).

## Scope of the emission model

UV1 = 200–400, UV2 = 250–400, UV3 = 300–400, VIS = 400–700 nm, with a constant
photon-counting throughput. Emission rates are interpolated in altitude and in
**the local SZA at each integration point**, Earth's shadow is tested at each
point, and the midpoint integral uses 2 km steps. A given satellite SZA does not
mean that the atmosphere along the line of sight has the same SZA. The band-mean
f_nu is computed from the photon radiance and converted to AB mag/arcsec².

The date sets the geometry; the emission rates remain the mean tables for
F10.7 = 100, Ap = 4, etc. The missing components of the tables, the absence of
UV absorption and of dawn/dusk chemistry history, and the global application of
the PALACE night layers still apply. The Moon direction includes the satellite
parallax; scattered moonlight is not included. For the GLOW terms of use, see
`docs/data/glowlicense.txt`.

## Zodiacal light + DGL + Airglow sum

Selecting **Zodiacal light + DGL + Airglow** in the web model list uses the same
orbit, phase, height, and ELV inputs as Airglow. In the Python GUI, choose
**Zodiacal + DGL + Airglow** under Map, or start in the summed mode:

```sh
python3 airglow_orbit.py --gui --model total
python3 airglow_orbit.py --model total --orbit noon-midnight \
  --orbital-phase-deg 90 --height 650 --date 2026-10-17T00:00:00Z \
  --band VIS --min-elv 5 --fits total_VIS_nside128_20261017.fits.gz \
  --png total_VIS_nside128_20261017.png
```

In the API use `sky.sky_map(..., model='total')`; `model='total'` is the default
and sums all three components. Airglow alone is `model='airglow'` or
`--model airglow` in the CLI. `--airglow-model combined/glow` selects the
airglow emission set independently of the summing mode.

The three components are **summed in MJy/sr in the same band** and then
converted to AB mag/arcsec². The photon radiance I_gamma is converted with
`I_airglow[MJy/sr] = 6.62607015e-10 * I_gamma / ln(lambda_max/lambda_min)`,
under the constant photon-counting throughput and the airglow AB definition.
Magnitudes are never added. The zodiacal light varies with the date (annual
harmonic templates of the SKYSURF model, combined at the Sun's J2000 ecliptic
longitude and rotated with the Sun), the DGL is fixed in RA/DEC, and the airglow
is integrated along the lines of sight of the selected orbit.

All three components use the same Sun/Moon avoidance, minimum ELV, and Earth
occultation. Directions where the airglow model is incomplete are also hatched
in the summed map. The FITS column INTENSITY_MJY_SR holds the total, and
ZODIACAL_MJY_SR, DGL_MJY_SR, AIRGLOW_MJY_SR, and AIRGLOW_PHOTON_RADIANCE hold
the components. The API `flux` is in MJy/sr in every mode. The airglow photon
unit is photons/s/cm²/sr, and `airglow_photon_radiance` keeps the airglow photon
radiance before summing.

## Display coordinates: Equatorial / Galactic

Select `equatorial` or `galactic` with the web **coordinate system** control or
the Python GUI **Coordinates** control. The default is equatorial RA, DEC; the
Galactic map shows longitude l and latitude b. The Sun, Moon, zenith, ELV and
avoidance boundaries, and contours are transformed to the same system. The web
page supports this for all five models.

```sh
python3 airglow_orbit.py --gui --coords galactic
python3 airglow_orbit.py --coords galactic --date 2026-10-17T00:00:00Z \
  --orbit noon-midnight --orbital-phase-deg 90 --model total \
  --fits total_VIS_nside128_20261017_galactic.fits.gz --png total_galactic.png
```

API: `sky.sky_map(date, coords='galactic', ...)`. The `sdrift_sky.py` CLI also
supports `--coords galactic` in every mode. Targets given with `--ra`, `--dec`
and `point()` queries are equatorial regardless of the map coordinate system.

Galactic pixel centers are converted to equatorial lines of sight, and the same
physical model is evaluated. The FITS file is therefore a RING HEALPix map in
the selected system, with `COORDSYS` C (equatorial) or G (Galactic). Masks and
components use the same pixel grid. Because the pixel centers differ between
frames, this is not a simple reordering of pixels.

The transformation uses the fixed matrix of the
[Astropy ICRS–Galactic rotation](https://docs.astropy.org/en/stable/coordinates/transforming.html),
identically in Python and on the web. It rotates the axes of the Sun and Moon
vectors and adds no parallax or aberration correction. The API `sun`, `moon`,
`state.zenith`, and RA/DEC metadata in the FITS header such as SUNRA/SUNDEC stay
in equatorial axes. Positions and pointer readouts in the GUI and the web page
use the selected system.

## Orbit options for Zodiacal / DGL / Zodiacal + DGL

On the web and in the Python orbit GUI, **every model** shares the orbit,
satellite height, phase, phase evolution from the reference UTC, minimum ELV,
and Sun/Moon avoidance settings. The defaults are H = 650 km, φ = 0°,
ELV ≥ 20°, Sun separation ≥ 90°, and Moon separation ≥ 60°. The zodiacal and
DGL brightness models are unchanged by the orbit; the orbit sets the observable
region and the zenith/ELV boundaries, and these modes neither integrate nor add
airglow. The airglow emission set applies only to the Airglow and total modes.

```sh
python3 airglow_orbit.py --gui --model dgl
python3 airglow_orbit.py --model zodiacal_dgl --coords galactic \
  --orbit noon-midnight --orbital-phase-deg 90 --height 650 --min-elv 20 \
  --date 2026-10-17T00:00:00Z --fits zodi_dgl_nside128_20261017_galactic.fits.gz
```

`--model` accepts `zodiacal`, `dgl`, `zodiacal_dgl`, `airglow`, and `total`.
The API equivalent is
`AirglowOrbitModel().sky_map(date, model='dgl', orbit='dawn-dusk', ...)`.
In these modes MODEL_VALID only indicates that the direction is not
Earth-occulted; it does not test airglow model validity or UV completeness.
Check `VIEW_ALLOWED` together with Earth occultation. `SkyModel` and both CLIs
support the same orbit options, and excluded directions are NaN in every mode.

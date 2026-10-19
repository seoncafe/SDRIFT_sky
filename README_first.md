# Building the S-DRIFT Dataset from Scratch

Last updated: 2026-10-20. The commands below are run from **the root of the
development folder that contains the generation scripts**. They build the
reference data for UV1, UV2, UV3 and VIS and convert them into the `docs/data/`
read by the web page and the common Python API. This document neither deletes
existing data nor starts any computation; run the commands for the steps you
need yourself.

**The published GitHub release alone cannot regenerate the data from the
original inputs.** The [SDRIFT_sky repository](https://github.com/seoncafe/SDRIFT_sky)
distributes the computed data and the viewing programs (`sdrift_sky.py`,
`airglow_orbit.py`, `airglow_orbit_gui.py`, `sky_coordinates.py`, `docs/`).
The generation scripts, tests, GLOW/PALACE/SKYSURF sources and the large
IRIS/SFD files named below are **in the development folder only**. To build
from scratch, prepare a backup of the full development folder and the original
inputs. To draw maps from the published data only, use the
[Python interface guide](README_python_interface.md).

## 1. Generation overview and bands

```text
SKYSURF source ── zodiacal_healpix.py ── zodi_data/*_template.fits.gz ─┐
IRIS + SFD + C2022_all ── dgl_healpix.py ── dgl_data/*_template.fits.gz ┤
GLOW ── prepare_airglow_glow.py ── glow_selected_lines_mean.npz ───────┤
                                    │                              │
PALACE ───────────────── prepare_airglow_palace.py                   │
                                    │                              │
                         glow_palace_night_mean.npz ────────────────┤
                                                                   ↓
                                                    export_web_data.py
                                                                   ↓
                                                    docs/data/manifest.json
                                                    docs/data/*.f32.gz
                                                                   ↓
                                          Common Python API / GUI / Web
```

Airglow is not computed by rotating a single all-sky reference map. The
emission rate `q(chi, h, component)` is stored first as a function of altitude
`h`, the solar zenith angle `chi` at the emitting point and the emission
component. Given an observation date, orbit, satellite altitude and orbital
phase, the table is interpolated and integrated along each line of sight to
produce a map.

| Band | Wavelength range |
|---|---|
| UV1 | 200–400 nm |
| UV2 | 250–400 nm |
| UV3 | 300–400 nm |
| VIS | 400–700 nm |

All bands assume constant throughput and photon-counting AB magnitudes.
The Zodiacal/DGL reference maps are RING HEALPix, NSIDE=128, equatorial
coordinates, with 196,608 pixels. The reference maps contain no Sun or Moon
avoidance mask. Avoidance conditions and the minimum ELV at the observation
time are applied when the final map is computed.

## 2. Development code and runtime environment

### 2.1 Required development files

Restore the following from a backup of the full development folder. The
original models alone are not enough to run the S-DRIFT generation and
conversion steps.

```text
SDRIFT_sky/
├── zodiacal_healpix.py
├── dgl_healpix.py
├── dgl_model/                    # modified DGL_lib.py, fitzpatrick.py, __init__.py
├── prepare_airglow_glow.py
├── prepare_airglow_palace.py
├── airglow_molecular.py
├── airglow_emission.py
├── airglow_geometry.py
├── export_web_data.py
├── sdrift_sky.py
├── airglow_orbit.py
├── airglow_orbit_gui.py
├── sky_coordinates.py
├── tests/                        # Python/JavaScript validation code
├── docs/assets/                  # web code
├── docs/index.html
├── skysurf/                      # original zodiacal light model (SKYSURF), Section 3
├── GLOW/                         # Fortran source including data/, Section 5
├── PALACE/                       # original archives and model/, Section 6
├── dgl_data/                     # original inputs of Section 4
├── zodi_data/                    # generated output
└── airglow_data/                 # generated output and build logs
```

Of this list, only `sdrift_sky.py`, `airglow_orbit.py`,
`airglow_orbit_gui.py`, `sky_coordinates.py` and `docs/` are part of the
published release. The optional H/SZA maps and LOS LUTs also need
`airglow_healpix.py`, `prepare_airglow_lut.py` and `airglow_lut.py`.
The validation scripts are listed in Section 10.

### 2.2 Python packages and compilers

Create a separate environment and use the same Python for installation and
generation. `python3` must refer to the Python of that environment.

```bash
python3 -m venv .venv-build
source .venv-build/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install "numpy>=2.0" scipy astropy healpy matplotlib Cython setuptools
mkdir -p zodi_data dgl_data airglow_data
```

Use Python 3.10 or newer. NumPy 2.0 or newer is required for generation (the
code uses `np.trapezoid()`). Cython is required even for the pure-Python
PALACE run. `zodiacal_healpix.py` imports the model folder directly, so no
separate zodisurf installation is needed. The published
`requirements-python.txt` covers only the Python interface and is not enough
for generation.

GLOW table generation needs an executable `gfortran` on PATH. Example
installation with macOS Homebrew:

```bash
brew install gcc
gfortran --version
```

Numba and MPI are not required. Everything can be generated with NumPy and
local processes. To use Numba, install a release compatible with the current
Python and NumPy and confirm that it imports. Do not assume that version
conflicts seen in an older environment apply to the current one.
`--engine numba` without Numba installed is an error for airglow/total maps.

```bash
python3 -m pip install numba
python3 -c "import numpy, numba; print(numpy.__version__, numba.__version__)"
```

The GUIs additionally need Tkinter/Tcl/Tk matching the Python version and a
graphical desktop. File generation steps do not need Tkinter.

## 3. Zodiacal light: annual-harmonic template FITS

### 3.1 Preparing the original model

The zodiacal light model is SKYSURF (O'Brien et al. 2025), a Kelsall-type
model, from `Zodiacal_Light_Model_Python/` in the
[SKYSURF repository](https://github.com/rosaliaobrien/skysurf). The local
directory name `skysurf/` is kept for compatibility.

Clone it only into a folder that does not have it yet; an existing copy need
not be cloned again. The commit used in the development folder is
`92dad09e18b9b69e51035a0fdf51364cb67e6ea6`.

```bash
git clone https://github.com/rosaliaobrien/skysurf.git skysurf
git -C skysurf checkout 92dad09e18b9b69e51035a0fdf51364cb67e6ea6
```

`skysurf/Zodiacal_Light_Model_Python/zodisurf/get_zmod.py` and `LICENSE` must be present.

### 3.2 Building the four bands

```bash
python3 zodiacal_healpix.py --make-template --nside 128 --band all \
  --reference-date 2000-03-20T07:35:00 --epochs 12 --workers 4 --no-plot
```

`--workers` is the number of local processes; adjust it to your resources.
`--epochs` is the number of model dates over one year used in the harmonic fit
(default 12, minimum 7). The build evaluates the model at each of these dates,
so it takes about 12 times as long as a single-date build: with NSIDE=128, all
bands and 12 workers it took roughly 2.5 h on the development machine.
This step writes all-sky FITS files without Sun or Moon masks.
`--no-plot` skips the optional PDF.

```text
zodi_data/zodiacal_UV1_nside128_template.fits.gz
zodi_data/zodiacal_UV2_nside128_template.fits.gz
zodi_data/zodiacal_UV3_nside128_template.fits.gz
zodi_data/zodiacal_VIS_nside128_template.fits.gz
```

SKYSURF brightness depends on the viewing direction relative to the Sun and,
through Earth's position relative to the tilted and offset interplanetary dust
cloud and the Earth–Sun distance, on the time of year. Each band template
therefore stores annual Fourier harmonics in the Sun's J2000 ecliptic
longitude `L`, sampled in a frame co-rotating with the Sun:

```text
I(v, t) = A0(u) + sum_{k=1..3} [Ck(u) cos(kL) + Sk(u) sin(kL)]
u = v rotated about the J2000 ecliptic pole by -(L - L_ref)
```

The FITS columns are `A0`, `C1`, `S1`, `C2`, `S2`, `C3`, `S3` (MJy/sr),
RING ordering, `COORDSYS='C'`. Header keywords: `REF-UTC` (frame epoch,
`2000-03-20T07:35:00`, the 2000 March equinox), `REFLON`, `HARMORD=3`,
`NEPOCH=12`, `FITRES` (maximum `|fit/model - 1|` over the fitted dates, all
sky), `FITRES3` (the same for Sun elongations of 3° or more),
`BANDID`/`BANDMIN`/`BANDMAX`, `MODEL='SKYSURF'` and `SOLARSP='LOGINTERP'`.
Older two-column templates are rejected with a message to rebuild them.

The 12 fitted dates are spaced evenly over one tropical year starting at
`REF-UTC`. Against direct model evaluations at dates within 2000 that were not
used in the fit, the maximum error is 0.33% and the median 0.01%; order 3 is
needed (order 2 leaves about 0.75%). At the fitted dates the maximum residual is 0.43% for Sun
elongations of 3° or more; within about 2° of the Sun it reaches 74–82%, far
inside the default 90° Sun avoidance. The annual phase uses the true J2000
ecliptic longitude of the Sun from Astropy. SKYSURF's internal solar-longitude
formula is tropical (of date) while its inputs are J2000 ecliptic, so
evaluating SKYSURF directly far from 2000 drifts by precession (about 1.4° per
century); the template, built from year-2000 dates, avoids that drift. Comparing the template with SKYSURF evaluated at the actual day therefore differs by up to about 0.04 mag at 1900 and 2100 (0.009 mag in 2026); against SKYSURF at the 2000 date with the same true solar longitude, the errors were measured on samples. The harmonic synthesis (UV1 and VIS, dates in 1900, 2026 and 2100, template pixels) stayed within about 0.003 mag. The full lookup in VIS, which adds rotation and NSIDE 128 bilinear interpolation, was evaluated on 2026 dates with up to 300 directions for each Sun-separation range: about 0.002 mag (largest 0.0021 mag) at 90 deg or more, the default avoidance; up to 0.004 mag at 10-90 deg; and up to about 0.06 mag at 3-10 deg, where interpolation of the steep brightness field dominates. These are sample measurements, not strict bounds for every date, band and direction. The regression tests in `tests/test_zodiacal_harmonics.py` use looser tolerances: 0.005 mag for the synthesis, and 0.005 / 0.005 / 0.01 / 0.08 mag for the lookup at >= 90 / 30-90 / 10-30 / 3-10 deg. The
remaining approximation is the perihelion drift relative to J2000 (about 0.3°
per century, < 0.1% in brightness).

The band spectral sampling step is 0.01 μm. The solar spectrum is interpolated
linearly in log(irradiance) between the nodes of SKYSURF's 0.05 μm table
instead of SKYSURF's nearest-node step function, and the model is evaluated
without SKYSURF's rounding to 1e-5 MJy/sr. Below 0.25 μm (the 200–250 nm part
of UV1) the SKYSURF albedo and phase-function fits are extrapolations of the
model calibration.

To build one band, specify e.g. `--band VIS`. The script overwrites existing
template files, so archive results you want to keep before running it.
The web converter reads the four NSIDE=128 file names above.

## 4. DGL: C2022_all reference FITS

### 4.1 Preparing the inputs

Place the following data in `dgl_data/`. They are all present in the
development folder; no symlink to `../DGL` or elsewhere is needed.

| File | Contents |
|---|---|
| `DGL_C2022_all.txt` | C2022 Fig. 3 right-panel correlation table; wavelength in Å |
| `IRIS_nohole_4_2048_v2.fits.gz` | IRIS 100 μm, Galactic RING, NSIDE=2048 |
| `lambda_sfd_ebv.fits.gz` | SFD E(B−V), Galactic NESTED, NSIDE=512 |
| `SOURCE_README.txt` | original description of the correlation table |
| `source_files.json` | source, copy path and SHA-256 of the inputs |

The correlation table and its description are also published as
`docs/data/DGL_C2022_all.txt` and `docs/data/DGL_source_README.txt`. IRIS and
SFD are available as HEALPix files (`IRIS_nohole_4_2048_v2.fits`,
`lambda_sfd_ebv.fits`) from
[NASA LAMBDA IRIS](https://lambda.gsfc.nasa.gov/product/foreground/fg_iris_get.html)
and [NASA LAMBDA SFD](https://lambda.gsfc.nasa.gov/product/foreground/fg_sfd_get.html);
gzip them after download. To reproduce the current results exactly, restore
**the same files from the backup**. Another processed IRIS map with the same
name, or a newly digitized table, is a different input.

Keep the project's modified copy in `dgl_model/`. Overwriting it with the
original `../DGL/DGL_lib.py` may change the wavelength-unit correction and the
self-absorption prescription that follows the paper. The code hashes in
`source_files.json` record the originals before copying and need not match the
current modified files. If you replace any data, update `source_files.json` to
match the actual inputs; its source fields are records, and the computation
never accesses those absolute paths.

### 4.2 Generation command and physical assumptions

```bash
python3 dgl_healpix.py --nside 128 --step-nm 0.5
```

Add `--overwrite` only when regenerating existing files.

```text
dgl_data/dgl_UV1_nside128_template.fits.gz
dgl_data/dgl_UV2_nside128_template.fits.gz
dgl_data/dgl_UV3_nside128_template.fits.gz
dgl_data/dgl_VIS_nside128_template.fits.gz
dgl_data/dgl_nside128_metadata.json
```

The output is RING HEALPix in equatorial coordinates with two columns,
`AB_MAG_ARCSEC2` and `INTENSITY_MJY_SR`. The original FIR and E(B−V) maps are
first degraded to NSIDE=128, the intensity is computed in Galactic coordinates
and interpolated to the equatorial pixel centers. The DGL itself does not
change with date.

The current prescription is:

```text
alpha_approx(lambda) = 2.1 × C2022_all(lambda)
beta = (1 − exp(−tau)) / tau
I_nu(lambda) = alpha_approx × beta × I100 × lambda_um / 100
<I_nu> = integral[I_nu(lambda) d(lambda)/lambda] / ln(lambda_max/lambda_min)
```

As in the self-absorption formula of
[Chellew et al. (2022)](https://arxiv.org/abs/2201.01378), **beta is
multiplied**. The original inverse-beta prescription is not used. Because the
table gives alpha_prime, the normalization approximation `alpha ≈ alpha_prime`,
`beta_ref=1` is used; the intrinsic alpha of the paper's sample is not exactly
recovered. Below 372.5 nm the correlation value at the shortest wavelength is
extended as a constant, while the Fitzpatrick extinction and the
wavelength-dependent intensity conversion are kept. The Fitzpatrick (FM_UNRED)
optical/IR part uses a natural cubic spline, as in the original
(`A(3000 Å)/E(B−V) = 5.5699`). The IRIS 100 μm map is used without subtracting
a zero level or the CIB; this is a deliberate choice of the project.

## 5. Airglow step 1: GLOW emission-rate table

### 5.1 Preparing the GLOW source

The original documentation is at [NCAR/HAO GLOW](https://www2.hao.ucar.edu/modeling/glow).
The generator is written for the Fortran file layout of the local
**GLOW v0.982**. To reproduce the results, restore this source and its `data/`
from the development folder. A newer GLOW release (for example from
[NCAR's GLOW repository](https://github.com/NCAR/GLOW)) must be checked for
compatibility with the build patches.

```text
GLOW/
├── Makefile.glowbasic.gfort
├── Glowlicense.txt
├── Glow.txt
├── cglow.f90, glow.f90, mzgrid.f90, gchem.f90, ...
└── data/
```

Running an existing `GLOW/glow.exe` is not needed.

### 5.2 Building the table

```bash
python3 prepare_airglow_glow.py --f107 100 --ap 4 --chi-step 1 \
  --top-km 1500 --radius-km 6371
```

Outputs:

```text
airglow_data/glow_selected_lines_mean.npz
airglow_data/glow_selected_lines_mean.json
airglow_data/build_glow/compile.log
airglow_data/build_glow/table.log
```

The main NPZ keys are `altitude_km`, `chi_deg`, `wavelength_nm`, `q` and
`metadata`. `q` has shape `(Nchi, Nz, Ncomponent)` in photons/cm³/s.
`q` is the volume photon production rate before division by 4π; the
line-of-sight integration applies `1/(4π)`. The components include individual
atomic lines and molecular channels treated at representative wavelengths.
The table version is `sdrift-glow-selected-lines-v3`.

The default conditions are F10.7=100 sfu, Ap=4, no auroral input and GLOW
xuvfac=3, with SZA 0–180° in 1° steps. The default GLOW altitude grid is
extended to 1500 km. For each SZA `chi`, the background atmosphere and
ionosphere (MSIS/IRI at the equator, longitude 0°) is the mean over four
seasons and the two local solar times `12 h ∓ chi/15 h` (the equinox, equator
relation; dawn and dusk branches). The night (chi=180°) columns are 630.0 nm
94.9 R, 636.4 nm 30.0 R, 557.7 nm 12.8 R and 520 nm 2.62 R. **These are
results for mean inputs, not a global climatological mean.**

The generator compiles patched sources in `airglow_data/build_glow/` and leaves
the originals unchanged. The patches are checked, and the build stops if a
patch does not apply. The `data` entry in this build is a relative symlink
(`../../GLOW/data`), so keep `GLOW/` and `airglow_data/` together in the
project root or rebuild. The default NPZ/JSON are overwritten on rerun.

## 6. Airglow step 2: combining PALACE nightglow components

### 6.1 Original archives and directories

The original distribution of PALACE v1.0 is
[Zenodo DOI 10.5281/zenodo.14064022](https://doi.org/10.5281/zenodo.14064022).
The development folder keeps the three archives below.

```text
PALACE/PALACE.zip
PALACE/test.zip
PALACE/PMD.zip
```

The model computation uses the code and data of `PALACE.zip`; `test.zip` holds
the distributed example for comparison. `PMD.zip` is auxiliary, but
`verify_palace.py` checks the CRC of all three archives, so **all three are
needed to run that verification.**

Unpack them only in a new environment where `model/` and `reference_test/` do
not exist yet, for example:

```python
from pathlib import Path
from zipfile import ZipFile

folder = Path("PALACE")
for archive, source, destination in [
    ("PALACE.zip", "PALACE", "model"),
    ("test.zip", "test", "reference_test"),
]:
    target = folder / destination
    if not target.exists():
        staging = folder / "unpacked"
        staging.mkdir(exist_ok=True)
        with ZipFile(folder / archive) as z:
            z.extractall(staging)
        (staging / source).rename(target)
```

The following files must then exist.

```text
PALACE/model/src/palace/palace.py
PALACE/model/src/palace/config/palace_default.par
PALACE/model/src/palace/data/palace_lines.fits
PALACE/model/src/palace/data/palace_cont.fits
PALACE/model/src/palace/data/palace_var.fits
PALACE/reference_test/palace_test.par
PALACE/reference_test/palace_test.fits
```

`prepare_airglow_palace.py` imports this local model directly. With Cython in
the virtual environment there is no need to build the extension modules with
`setup_cython.py`. If you bring an old Cython copy from `PALACE/dependencies/`
into the new environment, note that it takes import precedence over the
packages in the virtual environment.

### 6.2 Verifying the distributed example and combining

```bash
python3 verify_palace.py
python3 prepare_airglow_palace.py
```

If the combined table already exists, add `--overwrite` to the last command.
The verification writes `PALACE/verification.json` and
`airglow_data/palace_night_reference_300_700nm.fits.gz`. The reference spectrum
is a nighttime mean at the ground zenith, not a satellite all-sky map.

The combined outputs are:

```text
airglow_data/glow_palace_night_mean.npz
airglow_data/glow_palace_night_mean.json
airglow_data/glow_palace_night_mean_atomic_manifest.json
airglow_data/glow_palace_night_mean_molecular_manifest.json
```

The step adds PALACE OH/O2/Na, the continuum, non-overlapping O lines, the
mesospheric 557.7 nm residual and the 297.2 nm branch, and combines the
molecular bands from the same upper vibrational levels of N2/N2+.
The 557.7 nm residual is PALACE 162.67 R minus GLOW 12.83 R, i.e. 149.84 R
added as mesospheric 557.7 nm; the derived 297.2 nm branch is 6.36 R.
The molecular transition coefficients are in `airglow_molecular.py`, so no NIST
PDF is read at run time. PALACE column intensities are converted to assumed
Gaussian altitude profiles, and the additional nightglow components apply where
the emitting point lies in Earth's geometric shadow.

The default FWHMs are 8 km for OH and the green line, 10 km for other
mesospheric emissions and 100 km for thermospheric atomic lines. This is not a
full reconstruction of observed altitude profiles. PALACE below 300 nm, the H
geocorona, the full N2 LBH system, UV absorption and twilight chemical history
are not implemented. The latest component assessment is in
`AIRGLOW_COMPONENTS.md` (development folder only).

## 7. Optional: airglow FITS and LOS LUT

**Not required for the current common Python API or the web viewer.**
The web converter in Section 8 reads the two emission-rate NPZ files, not the
maps or `*_los_lut.npz` below. The current orbit viewer interpolates the band
emission-rate table and integrates along each line of sight.

### 7.1 Reference maps with fixed H/SZA

Use these to check results in the older reference geometry or to archive them
as FITS.

```bash
python3 airglow_healpix.py --table airglow_data/glow_selected_lines_mean.npz \
  --H 550 --sza 90 --band UV1 UV2 UV3 VIS --nside 128 --workers 4 --pdf-dir output/pdf/airglow --contour-step 0.2
python3 airglow_healpix.py --table airglow_data/glow_palace_night_mean.npz \
  --H 550 --sza 90 --band UV1 UV2 UV3 VIS --nside 128 --workers 4 --pdf-dir output/pdf/airglow --contour-step 0.2
```

For example, the VIS outputs are:

```text
airglow_data/airglow_VIS_nside128_template_H550_sza90.fits.gz
airglow_data/airglow_VIS_nside128_template_H550_sza90_glow_palace.fits.gz
```

These maps use a reference geometry with fixed H and SZA and do not represent
an orbital position on a specific date. H=550 km above is a test case; the
current orbit viewer default altitude is 650 km. The `airglow_healpix.py` CLI
contour default is 0.5 mag, so the commands set 0.2 explicitly. The FITS `HMIN`
is the lowest supported altitude of the table (80 km for the combined table).
Below that altitude results are flagged in `MODEL_VALID`; above the table top
(1500 km) emission is taken as zero without a flag. Add `--overwrite` to
recompute existing files.

### 7.2 Integrated LUT for the older H/SZA viewer

```bash
python3 prepare_airglow_lut.py --table airglow_data/glow_selected_lines_mean.npz --workers 4
python3 prepare_airglow_lut.py --table airglow_data/glow_palace_night_mean.npz --workers 4
```

These write `glow_selected_lines_mean_los_lut.npz` and
`glow_palace_night_mean_los_lut.npz` in the same folder.
The LUT integrates many H/SZA/azimuth/tangent-altitude combinations, so this
step can take a long time. Rebuild the LUT whenever its emission-rate table
changes; an existing LUT cannot be used if the source SHA-256 differs. Add
`--overwrite` to regenerate.

Compare the LUT error with direct integration:

```bash
python3 validate_airglow_lut.py --table airglow_data/glow_palace_night_mean.npz --nside 16
```

You can go straight to Section 8; the common Python API and web maps work
without these files.

## 8. Converting to shared Python/web data

Run this after all of the following are ready.

- `zodi_data/zodiacal_{UV1,UV2,UV3,VIS}_nside128_template.fits.gz`
- `dgl_data/dgl_{UV1,UV2,UV3,VIS}_nside128_template.fits.gz`
- `dgl_data/dgl_nside128_metadata.json`, the correlation table and its description
- `airglow_data/glow_selected_lines_mean.npz`
- `airglow_data/glow_palace_night_mean.npz`
- `GLOW/Glowlicense.txt`, `skysurf/Zodiacal_Light_Model_Python/LICENSE`

```bash
python3 export_web_data.py
```

This step generates:

| Output | Purpose |
|---|---|
| `docs/data/manifest.json` | bands, dimensions, grids, sources, SHA-256, model assumptions |
| `zodiacal_{band}.f32.gz` | Zodiacal harmonic coefficients, shape `[7, 196608]`, MJy/sr |
| `dgl_{band}.f32.gz` | DGL reference map in MJy/sr |
| `airglow_{glow,combined}_{band}.f32.gz` | `(Nchi, Nz, 2)` band emission-rate table |
| `sun_j2000.f32.gz` | Sun direction every 12 hours |
| `moon_gcrs.f32.gz` | geocentric GCRS Moon vectors in km every 3 hours |
| license and source files | GLOW, Zodiacal and DGL data descriptions |

All file names in the table are under `docs/data/`, and `{band}` stands for
each of the four bands. The files are gzipped little-endian Float32 arrays;
`.f32.gz` is not FITS. The manifest entry `zodiacal[band]` carries
`harmonic_order` and `harmonics` (`A0` … `S3`); Python and the browser combine
the harmonics for the requested date and then sample the rotated map. The Moon
table has `lunar.units = 'km'`, so the orbit models can apply the parallax of
the satellite position. The two values on the last airglow axis are the band
photon production rates of the regular component and of the nightglow
component applied only inside Earth's shadow; they are not yet
line-of-sight-integrated radiances. The Sun and Moon tables are generated with
the Astropy builtin ephemeris for 1900–2100. The full time range is computed,
so the conversion can take a while.

**The exporter can finish without DGL if the DGL inputs are missing.**
To publish the full model, generate DGL first and confirm in Section 10 that
the manifest contains `dgl`. The converter uses fixed NSIDE=128 file names, so
files built only at another resolution cannot be converted with this command.

If only DGL changed, the existing manifest can be updated in place:

```bash
python3 export_web_data.py --dgl-only
```

This option and `--moon-only` require an existing manifest. For the first run
from the original inputs, use the full conversion without options. Use the full
conversion also after Zodiacal or airglow changes. The exporter overwrites its
outputs and the manifest, so generate and validate in the development folder
before deploying.

## 9. Final maps and summed models

Use the current common interface (`sdrift_sky.py` is also in the published
release). To avoid confusion with older CLI options of the reference-data
scripts, the final observation map is selected with `--model` and the airglow
component model with `--airglow-model`.

```bash
python3 sdrift_sky.py --model total --airglow-model combined --band VIS \
  --date 2026-10-19T00:00:00Z --orbit dawn-dusk --height 650 \
  --orbital-phase-deg 0 --min-elv 20 --sun-separation 90 --moon-separation 60 \
  --nside 128 --contour-step 0.2 \
  --fits output/fits/maps/total_VIS_nside128_20261019.fits.gz \
  --pdf output/pdf/maps/total_VIS_nside128_20261019.pdf
```

`total` is Zodiacal + DGL + airglow. The three components are **added in MJy/sr
in the same band** and then converted to AB mag/arcsec²; magnitudes are never
added. The final output uses the 14 FITS columns of common API v2. Do not
confuse it with the reference template FITS columns or the older H/SZA map
format.

The orbit is a circular sun-synchronous orbit of radius R+H (R = 6371 km). The
inclination follows from the J2 sun-synchronous condition (97.96° at H = 650 km),
and the ascending node is at a fixed mean local solar time: LTAN 18:00 for
dawn-dusk and 12:00 for noon-midnight. The orbital phase follows the
`NORTH_V2` convention: 0° is the northernmost point, 90° the descending node
(dawn side for dawn-dusk, midnight side for noon-midnight), 180° the
southernmost point and 270° the ascending node; the phase advances along the
motion. In epoch mode (`--epoch`) the phase advances with the draconic period
(97.70 min at 650 km). The satellite SZA varies around the orbit; near the June
solstice a dawn-dusk orbit at 650 km enters Earth's shadow. The Moon direction
is the one seen from the satellite. FITS headers record `ORBMODEL`, `INCL`,
`RAAN`, `LTAN`, `PHCONV='NORTH_V2'`, `PERIOD`, `MOONOBS='SATELLITE'`, and
`EPOCH` only in epoch mode. See
[README_airglow_orbit.md](README_airglow_orbit.md) for details.

Example of querying selected coordinates. RA/DEC are in degrees and the date
is UTC. Outside the observing conditions the brightness is NaN, so check
`visible` as well.

```python
from sdrift_sky import SkyModel

sky = SkyModel()
p = sky.point(
    30, 10, "2026-10-19T00:00:00Z", model="total", band="VIS",
    orbit="noon-midnight", orbital_phase_deg=90,
)
print(p["ab_mag_arcsec2"], p["intensity_mjy_sr"], p["visible"])
```

Checking the web view:

```bash
python3 -m http.server 8765 --directory docs
```

Open `http://localhost:8765/` in a browser. Restart any running Python/GUI
session with the new model and reload the browser so that cached old data are
not mixed in.

## 10. Validating the results

### 10.1 Recording sources, outputs and environment

After generation, keep the original archives, code version, Python package
versions, command lines and logs together with the JSON metadata of the
FITS/NPZ files. A different SHA-256 does not always mean a different physical
result, but it must not be assumed to be the same input as the previous data.

```bash
python3 -m pip freeze > airglow_data/build_python_packages.txt
```

The following example checks that the generated files exist, their FITS
resolution and coordinate system, and the exporter hashes. Note that the hashes
are of the uncompressed arrays.

```python
from pathlib import Path
import gzip
import hashlib
import json
import numpy as np
from astropy.io import fits

root = Path(".")
zodi_columns = ["A0", "C1", "S1", "C2", "S2", "C3", "S3"]
for band in ["UV1", "UV2", "UV3", "VIS"]:
    for directory, prefix, columns in [
        ("zodi_data", "zodiacal", zodi_columns),
        ("dgl_data", "dgl", ["INTENSITY_MJY_SR"]),
    ]:
        path = root / directory / f"{prefix}_{band}_nside128_template.fits.gz"
        with fits.open(path, checksum=True) as hd:
            hd.verify("exception")
            header = hd[1].header
            assert header["NSIDE"] == 128
            assert header["ORDERING"] == "RING"
            assert header["COORDSYS"] == "C"
            if prefix == "zodiacal":
                assert header["HARMORD"] == 3 and header["NEPOCH"] >= 7
            for name in columns:
                values = hd[1].data[name].ravel()
                assert values.size == 12 * 128**2
                assert np.all(np.isfinite(values))
            if prefix == "dgl":
                assert np.all(values >= 0)

manifest = json.loads((root / "docs/data/manifest.json").read_text())
assert set(manifest["bands"]) == {"UV1", "UV2", "UV3", "VIS"}
assert "dgl" in manifest
assert set(manifest["airglow"]) == {"glow", "combined"}
assert manifest["lunar"]["units"] == "km"

def verify_arrays(item):
    if isinstance(item, dict):
        if "file" in item and "sha256" in item:
            raw = gzip.decompress((root / "docs/data" / item["file"]).read_bytes())
            assert hashlib.sha256(raw).hexdigest() == item["sha256"], item["file"]
            assert len(raw) == int(np.prod(item["shape"])) * 4, item["file"]
        for value in item.values():
            verify_arrays(value)
    elif isinstance(item, list):
        for value in item:
            verify_arrays(value)

verify_arrays(manifest)
print("PASS: templates, all model groups, exported array sizes and SHA-256")
```

If a FITS file has no CHECKSUM, `checksum=True` alone does not guarantee one is
present. The code above also checks the FITS structure and values, and the
converted arrays are checked with SHA-256.

### 10.2 Local validation scripts

The tests live only in the development folder, in the `tests/` package (see
`tests/README.md` there); they are not part of the published release. Run them
from the root of the development folder. Tests that compare with the
JavaScript implementation call Node.js from Python, so Node.js is required.
The tests need the generated `docs/data/`, the reference FITS and the airglow
emission-rate NPZ files.

```bash
python3 -m unittest tests.test_dgl tests.test_airglow tests.test_airglow_components
python3 -m tests.test_web_sky
python3 -m unittest tests.test_airglow_orbit tests.test_sky_coordinates tests.test_python_interface tests.test_zodiacal_harmonics tests.test_orbit_gui_state
python3 -m tests.test_moon_avoidance
node tests/test_web_dgl.mjs
node tests/test_web_moon.mjs
node tests/test_web_dates.mjs
node tests/test_web_worker_requests.mjs
node tests/test_web_app_cancel.mjs
```

The Python/JavaScript parity tests also cover the sun-synchronous orbit
elements, the Moon parallax and the zodiacal harmonics.

Only if the LOS LUTs were built, also run:

```bash
python3 -m unittest tests.test_airglow_lut
```

These tests check numerical integration, units, masks and agreement between
implementations; they do not guarantee the completeness of all components or
the physical accuracy relative to real satellite observations.

## 11. What to regenerate after changing inputs

| Changed | Regenerate |
|---|---|
| Zodiacal source, reference UTC, number of fit epochs, spectral sampling | Zodiacal FITS → full export |
| IRIS/SFD, correlation table, DGL prescription | DGL FITS/metadata → DGL-only or full export |
| GLOW source, F10.7, Ap, SZA step, altitude top | GLOW NPZ → combined PALACE NPZ → full export |
| PALACE source, layer widths, component selection, molecular branches | combined NPZ → full export |
| LUTs for the older H/SZA viewer after the airglow changes above | those LUTs and map caches |
| Observation UTC, orbit, phase, H, ELV, avoidance angles | no source regeneration; recompute the final map |

If `docs/data/`, which the common API reads, is not updated, rebuilding the NPZ
or reference FITS does not change the displayed values. Update the manifest and
its arrays together.

## 12. Files to publish on GitHub

Users need the updated `docs/data/`, the web code, the common Python code, the
requirements file and the guides. With the development folder and the public
clone kept separate, copy `docs/data/` to `web_publish/docs/data/` after
validation and publish from there. Adding or editing this document alone does
not require regenerating or recopying the data.

The original GLOW/PALACE sources, the large IRIS/SFD FITS files, build
products, older maps and logs are not needed by users and are not included in
the user release. Keep `.gitignore` out of the pushed files, as before. Back up
the full development environment for reproduction separately from the public
release.

Preserve the license and attribution of each source. Use of GLOW calculations
is governed by the GLOW academic research license; the PALACE code is GPLv3
and its distributed model data are CC BY 4.0. This project is software that
displays results computed with GLOW and the other models.

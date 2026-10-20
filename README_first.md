# Rebuilding the S-DRIFT Dataset from Scratch

Last updated: 2026-10-20.

This document is the **handover procedure for a new user who has received the
full development folder and wants to regenerate the dataset for every band
from the original inputs**. It does not rely on computed results already in
the folder. It explains how to build the Zodiacal/DGL reference FITS files,
the GLOW/PALACE emission-rate tables, the optional LOS LUTs, the final maps
and the data shared by the web and Python interfaces. Section 9.1 also explains
how to compute sensitivity maps from the three background components, using
the default f_sys = 0.1%, usable exposure = 10 hours and 5×5 pixel binning.

Copy the received folder to a new working location and run everything from its
root. The current folder already contains all original model code and DGL
inputs. If you received it intact, use the existing files. If something is
missing, prepare it following the download or backup-restore instructions in
Sections 3.1, 4.1, 5.1 and 6.1, then rerun the input check.
With all original inputs present, the procedure does not depend on another
user's home directory, `../DGL`, an existing `docs/data/` or additional GitHub
clones. Existing inputs need not be downloaded again. Internet access is needed
only for missing inputs, Python packages and compilers. No GitHub commit/push is
performed.

The order is **Section 2 environment and input check → Sections 3–6 source
computations → Section 7 optional computations → Section 8 conversion to shared
data → Sections 9–10 running and validation**.
Section 13 contains a single command block to run after the environment is ready.
To query existing data only, see the
[Python interface guide](README_python_interface.md). The commands here are
instructions; large recomputations were not run while writing this document.

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

The full folder to hand over contains the code and source data below.
The recipient keeps this directory structure and works from its root.

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
├── sensitivity.py                # CCD detection limit from the summed sky background
├── airglow_orbit.py
├── airglow_orbit_gui.py
├── sky_coordinates.py
├── figs/                         # geometry figure script and dawn-dusk PDF build script
├── output/tex/                   # TeX sources of the dawn-dusk documents
├── tests/                        # Python/JavaScript validation code
├── docs/assets/                  # web code
├── docs/index.html
├── skysurf/                      # original zodiacal light model (SKYSURF)
├── GLOW/                         # Fortran source including data/
├── PALACE/                       # original archives and model/
├── dgl_data/                     # original FITS, correlation table, source JSON
├── zodi_data/                    # generated output
└── airglow_data/                 # generated output and build logs
```

A handover for full regeneration also includes `airglow_healpix.py`,
`prepare_airglow_lut.py`, `airglow_lut.py`, `validate_airglow_lut.py`,
`verify_palace.py`, the tests of Section 10 and the `tests/test_*.mjs` files
they call. Include `README_first.md`, `README_python_interface.md`, `README_sensitivity.md`
and the requirements files as well. The sensitivity web calculation also needs
`docs/assets/sensitivity.mjs`; include `tests/test_sensitivity.py` and
`tests/test_sensitivity.mjs` for validation.

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
separate zodisurf installation is needed. To rebuild
`output/pdf/airglow_geometry.pdf` (Section 14) also install `reportlab`.

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

### 2.3 Checking the received sources and tools

If inputs are missing, first apply the recovery instructions in Sections 3.1,
4.1, 5.1 and 6.1. In a terminal, move to the project root and run the block
below. It checks for missing input files and the hashes of the DGL source data,
and does not refer to paths on the previous computer. Install `gfortran` and
`node` first if they are absent. Node.js is used for the web-implementation
comparison tests and is not needed for the Python computations themselves.

```bash
gfortran --version
node --version
python3 - <<'PY'
from pathlib import Path
import hashlib
import json
import shutil

root = Path.cwd()
required = [
    'zodiacal_healpix.py', 'dgl_healpix.py', 'dgl_model/__init__.py',
    'dgl_model/DGL_lib.py', 'dgl_model/fitzpatrick.py',
    'prepare_airglow_glow.py', 'prepare_airglow_palace.py',
    'airglow_molecular.py', 'airglow_emission.py', 'airglow_geometry.py',
    'airglow_healpix.py', 'airglow_lut.py', 'prepare_airglow_lut.py',
    'validate_airglow_lut.py', 'verify_palace.py', 'export_web_data.py',
    'sdrift_sky.py', 'airglow_orbit.py', 'airglow_orbit_gui.py', 'sky_coordinates.py',
    'sensitivity.py', 'docs/assets/sensitivity.mjs', 'README_sensitivity.md',
    'tests/test_sensitivity.py', 'tests/test_sensitivity.mjs',
    'skysurf/Zodiacal_Light_Model_Python/zodisurf/get_zmod.py',
    'skysurf/Zodiacal_Light_Model_Python/LICENSE',
    'GLOW/Makefile.glowbasic.gfort', 'GLOW/Glowlicense.txt', 'GLOW/data',
    'GLOW/cglow.f90', 'GLOW/glow.f90', 'GLOW/mzgrid.f90',
    'PALACE/PALACE.zip', 'PALACE/test.zip', 'PALACE/PMD.zip',
    'dgl_data/DGL_C2022_all.txt', 'dgl_data/IRIS_nohole_4_2048_v2.fits.gz',
    'dgl_data/lambda_sfd_ebv.fits.gz', 'dgl_data/source_files.json',
    'dgl_data/SOURCE_README.txt', 'docs/index.html', 'docs/assets',
    'figs/draw_airglow_geometry.py', 'figs/build_dawn_dusk_document.sh',
    'output/tex/dawn_dusk_geometry.tex', 'output/tex/dawn_dusk_airglow_geometry.tex',
    'tests/__init__.py', 'tests/test_dgl.py', 'tests/test_airglow.py', 'tests/test_airglow_components.py',
    'tests/test_airglow_lut.py', 'tests/test_web_sky.py', 'tests/test_web_sky.mjs',
    'tests/test_airglow_orbit.py', 'tests/test_airglow_orbit.mjs',
    'tests/test_sky_coordinates.py', 'tests/test_python_interface.py',
    'tests/test_moon_avoidance.py', 'tests/test_web_dgl.mjs', 'tests/test_web_moon.mjs',
    'tests/test_zodiacal_harmonics.py', 'tests/test_zodiacal_harmonics.mjs',
    'tests/test_web_dates.mjs', 'tests/test_web_worker_requests.mjs',
    'tests/test_web_app_cancel.mjs', 'tests/test_orbit_gui_state.py',
]
missing = [name for name in required if not (root / name).exists()]
if missing:
    raise SystemExit('Missing from handover:\n' + '\n'.join(missing))
for tool in ('gfortran', 'node'):
    if shutil.which(tool) is None:
        raise SystemExit(f'{tool} not found on PATH.')
for item in json.loads((root / 'dgl_data/source_files.json').read_text()):
    relative = item['copy']
    # Source-code hashes may differ from the modified copies; check input data only.
    if relative.startswith('dgl_data/'):
        actual = hashlib.sha256((root / relative).read_bytes()).hexdigest()
        if actual != item['sha256']:
            raise SystemExit(f'DGL input hash mismatch: {relative}')
print('PASS: supplied inputs, generation/test code, gfortran and Node.js')
PY
```

This check covers representative files and the DGL input hashes. The full GLOW
source with `data/` and the SKYSURF model and auxiliary data must be handed
over as whole folders. Include the actual input files rather than symlinks to
external folders.

### 2.4 Preserving old results and starting fresh

Run the block below once in a **working copy, not in the received original
folder**. It moves existing results and the PALACE modules of the previous
environment to `rebuild_archive/<UTC>/` and creates empty output directories.
The input DGL FITS files, the correlation table and the model archives are
kept. Old symlinks and binaries in the GLOW build folder are moved as well.

```bash
python3 - <<'PY'
from datetime import datetime, timezone
from pathlib import Path
import shutil

root = Path.cwd()
backup = root / 'rebuild_archive' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
backup.mkdir(parents=True, exist_ok=False)

def preserve(relative):
    source = root / relative
    if source.exists() or source.is_symlink():
        destination = backup / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(source), str(destination))

for relative in ['zodi_data', 'airglow_data', 'docs/data',
                 'PALACE/model', 'PALACE/reference_test', 'PALACE/dependencies']:
    preserve(relative)
for pattern in ['dgl_*_nside*_template.fits.gz', 'dgl_nside*_metadata.json']:
    for source in sorted((root / 'dgl_data').glob(pattern)):
        preserve(source.relative_to(root))
for relative in ['zodi_data', 'airglow_data', 'docs/data', 'rebuild_logs',
                 'output/pdf/zodiacal', 'output/pdf/airglow',
                 'output/pdf/maps', 'output/fits/maps']:
    (root / relative).mkdir(parents=True, exist_ok=True)
print('Previous files archived in:', backup)
PY
```

The PALACE archives are then unpacked again in Section 6.1; that code is also
part of the Section 13 block. Do not rerun the cleanup block in the middle of
a rebuild. If the Section 13 run fails, fix the cause and resume from the
failed step. Apply `--overwrite` only to the outputs needed for resuming.

Record package and compiler versions as well.

```bash
python3 -m pip freeze > rebuild_logs/python_packages.txt
python3 --version > rebuild_logs/python_version.txt
gfortran --version > rebuild_logs/gfortran_version.txt
node --version > rebuild_logs/node_version.txt
```

## 3. Zodiacal light: annual-harmonic template FITS

### 3.1 Preparing the original model

The zodiacal light model is SKYSURF (O'Brien et al. 2025), a Kelsall-type
model, from `Zodiacal_Light_Model_Python/` in the
[SKYSURF repository](https://github.com/rosaliaobrien/skysurf). The local
directory name `skysurf/` is kept for compatibility.

`skysurf/` already contains the code and model data, so use it as is.
**If the folder is missing**, download it from the original repository as below.
If only some files are missing from an existing folder, restore them by
comparison with a fresh copy in a separate folder.

```bash
# Only if skysurf/ is missing
git clone https://github.com/rosaliaobrien/skysurf.git skysurf
git -C skysurf checkout 92dad09e18b9b69e51035a0fdf51364cb67e6ea6
```

The original commit of the supplied copy is
`92dad09e18b9b69e51035a0fdf51364cb67e6ea6`. A missing `.git/` does not affect
the model computation.
`skysurf/Zodiacal_Light_Model_Python/zodisurf/get_zmod.py` and `LICENSE` must be present.

### 3.2 Building the four bands

```bash
python3 zodiacal_healpix.py --make-template --nside 128 --band all \
  --reference-date 2000-03-20T07:35:00 --epochs 12 --workers 4 --no-plot
```

`--workers` is the number of local processes; adjust it to your resources.
`--epochs` is the number of model dates over one year used in the harmonic fit
(default 12, minimum 7). The template build evaluates the model at each of these
dates, so it takes about 12 times as long as a single-date build: with NSIDE=128,
all bands and 12 workers it took roughly 2.5 h on the development machine.
This step writes all-sky FITS files without Sun or Moon masks.
`--no-plot` skips the optional PDF.

```text
zodi_data/zodiacal_UV1_nside128_template.fits.gz
zodi_data/zodiacal_UV2_nside128_template.fits.gz
zodi_data/zodiacal_UV3_nside128_template.fits.gz
zodi_data/zodiacal_VIS_nside128_template.fits.gz
```

SKYSURF brightness depends on the viewing direction relative to the Sun and,
through Earth's position relative to the tilted and offset interplanetary
dust cloud and the Earth–Sun distance, on the time of year. Each band template
therefore stores annual Fourier harmonics in the Sun's J2000 ecliptic
longitude `L`, sampled in a frame co-rotating with the Sun:

```text
I(v, t) = A0(u) + sum_{k=1..3} [Ck(u) cos(kL) + Sk(u) sin(kL)]
u = v rotated about the J2000 ecliptic pole by -(L - L_ref)
```

The FITS columns are `A0`, `C1`, `S1`, `C2`, `S2`, `C3`, `S3` (MJy/sr),
RING ordering, `COORDSYS='C'`. Header keywords: `REF-UTC` (frame epoch,
`2000-03-20T07:35:00`, the 2000 March equinox), `REFLON` (Sun J2000 ecliptic
longitude at `REF-UTC`), `HARMORD=3`, `NEPOCH=12`, `FITRES` (maximum
`|fit/model - 1|` over the fitted dates, all sky), `FITRES3` (the same for
Sun elongations of 3° or more), `BANDID`/`BANDMIN`/`BANDMAX`,
`MODEL='SKYSURF'` and `SOLARSP='LOGINTERP'`. Older two-column templates are
rejected with a message to rebuild them.

The 12 fitted dates are spaced evenly over one tropical year starting at
`REF-UTC`; the seven coefficient maps are least-squares fits. Against direct
model evaluations at dates within 2000 that were not used in the fit, the
maximum error is 0.33% and the median 0.01%; order 3 is needed (order 2 leaves
about 0.75%). For Sun elongations of 3° or more the maximum residual at the fitted dates is 0.43% (header FITRES3; for the current templates from a 2000-pixel sample). Within about 2° of the Sun it reaches 74–82% (FITRES, all sky), because SKYSURF's internal solar longitude differs from the true one by up to 0.024° where the brightness gradient is extreme; that region lies far inside the default 90° Sun avoidance.

The annual phase uses the true J2000 ecliptic longitude of the Sun from
Astropy. SKYSURF's internal solar-longitude formula is tropical (of date)
while its longitude/latitude inputs are J2000 ecliptic, so evaluating SKYSURF
directly far from 2000 drifts by precession (about 1.4° per century). The
template, built from year-2000 dates where both agree to better than 0.03°,
avoids that drift. Comparing the template with SKYSURF evaluated at the actual day therefore differs by up to about 0.04 mag at 1900 and 2100 (0.009 mag in 2026); against SKYSURF at the 2000 date with the same true solar longitude, the errors were measured on samples. The harmonic synthesis (UV1 and VIS, dates in 1900, 2026 and 2100, template pixels) stayed within about 0.003 mag. The full lookup in VIS, which adds rotation and NSIDE 128 bilinear interpolation, was evaluated on 2026 dates with up to 300 directions for each Sun-separation range: about 0.002 mag (largest 0.0021 mag) at 90 deg or more, the default avoidance; up to 0.004 mag at 10-90 deg; and up to about 0.06 mag at 3-10 deg, where interpolation of the steep brightness field dominates. These are sample measurements, not strict bounds for every date, band and direction. The regression tests in `tests/test_zodiacal_harmonics.py` use looser tolerances: 0.005 mag for the synthesis, and 0.005 / 0.005 / 0.01 / 0.08 mag for the lookup at >= 90 / 30-90 / 10-30 / 3-10 deg. The remaining approximation is the perihelion drift
relative to J2000 (about 0.3° per century, < 0.1% in brightness).

Other details:

- The band spectral sampling step is 0.01 μm.
- SKYSURF's `solar_sp` returns the nearest node of a 0.05 μm table, a step
  function with jumps up to ×10 between 0.20 and 0.25 μm. The project keeps the
  node values and interpolates log(irradiance) linearly between nodes
  (`zodiacal_healpix.solar_irradiance`). Band means with the old step function
  were biased high by 6.4% (UV1), 5.7% (UV2), 4.0% (UV3) and 0.3% (VIS).
- SKYSURF's `get_zmod` rounds its result to 1e-5 MJy/sr, comparable to the UV
  signal; `zodiacal_healpix.zodiacal_model` calls the same kernel and
  parameters without rounding.
- Below 0.25 μm (the 200–250 nm part of UV1) the SKYSURF albedo and
  phase-function fits are extrapolations of the model calibration.

To build one band, specify e.g. `--band VIS`. The script overwrites existing
template files, so archive results you want to keep before running it.
The web converter reads the four NSIDE=128 file names above.

## 4. DGL: C2022_all reference FITS

### 4.1 Preparing the inputs

Place the following data in `dgl_data/`. They are all already copied into the
current development folder. No symlink to `../DGL` or elsewhere is needed.

| File | Contents |
|---|---|
| `DGL_C2022_all.txt` | C2022 Fig. 3 right-panel correlation table; wavelength in Å |
| `IRIS_nohole_4_2048_v2.fits.gz` | IRIS 100 μm, Galactic RING, NSIDE=2048 |
| `lambda_sfd_ebv.fits.gz` | SFD E(B−V), Galactic NESTED, NSIDE=512 |
| `SOURCE_README.txt` | original description of the correlation table |
| `source_files.json` | source, copy path and SHA-256 of the inputs |

These five files and the `dgl_model/` code are already in the folder. Do not
download them again if they are present. **If any are missing**, download them
from the sources below or restore them from a project backup under the file
names above.

| Missing item | Download or restore |
|---|---|
| IRIS 100 μm | From [NASA LAMBDA IRIS](https://lambda.gsfc.nasa.gov/product/foreground/fg_iris_get.html), choose HEALPix `IRIS_nohole_4_2048_v2.fits`; NSIDE=2048 |
| SFD E(B−V) | From [NASA LAMBDA SFD](https://lambda.gsfc.nasa.gov/product/foreground/fg_sfd_get.html), choose HEALPix `lambda_sfd_ebv.fits`; NSIDE=512 |
| C2022_all table | Save the [correlation table kept in S-DRIFT](https://github.com/seoncafe/SDRIFT_sky/blob/main/docs/data/DGL_C2022_all.txt) as a raw file or restore it from the full handover backup |
| Original table description | Save the [description kept in S-DRIFT](https://github.com/seoncafe/SDRIFT_sky/blob/main/docs/data/DGL_source_README.txt) as `SOURCE_README.txt` or restore it from a backup |
| Modified `dgl_model/`, `source_files.json` | Restore from the full development folder backup or the provider; no separate public download address is known |

For IRIS/SFD, download the **HEALPix files**, not the Mollweide or ZEA image
files. Put the uncompressed FITS in `dgl_data/` and gzip only the needed files.
Do not merely rename the extension to `.gz`.

```bash
# Run each only if the .fits.gz is absent and the downloaded .fits is present
gzip -n dgl_data/IRIS_nohole_4_2048_v2.fits
gzip -n dgl_data/lambda_sfd_ebv.fits
```

Check the resolution, ordering, units and source of downloaded data. A
different gzip method or header can change the SHA-256 relative to the existing
compressed file, so after confirming the input is the same, record the new
download URL, version and hash in `source_files.json`, keeping the original
record. For a bit-identical reproduction of the supplied files, restore the
same files from a backup.

Absolute paths of the previous computer in the `source` field of
`source_files.json` are provenance records; the computation does not use them.
To reproduce the current results, use the supplied input files.
Another processed IRIS map with the same name, or a newly digitized table, is a
different input.

Keep the project's modified copy in `dgl_model/`. Overwriting it with the
original `../DGL/DGL_lib.py` may change the wavelength-unit correction and the
self-absorption prescription that follows the paper. The code hashes in
`source_files.json` record the originals before copying and need not match the
current modified files.

If you replace any data, update `source_files.json` to match the actual inputs.
The source fields are records; the computation never accesses those absolute paths.

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
**GLOW v0.982**. `GLOW/` already contains the source and the input `data/` of
this version, so use it as is. **If the folder is missing**, download the
source together with `data/` from [NCAR's GLOW repository](https://github.com/NCAR/GLOW).

```bash
# Only if GLOW/ is missing
git clone https://github.com/NCAR/GLOW.git GLOW
```

This clone fetches the repository's default version, which is not guaranteed
to be the supplied v0.982. Check `Glow.txt`, the required files and build
compatibility, and record the code version. To reproduce the supplied results
exactly, restore the same v0.982 source and data from a backup. Do not
silently replace it with the latest version.

```text
GLOW/
├── Makefile.glowbasic.gfort
├── Glowlicense.txt
├── Glow.txt
├── cglow.f90, glow.f90, mzgrid.f90, gchem.f90, ...
└── data/
```

The official page and the NCAR repository are the sources for obtaining the
code and checking its license. If the current folder is complete, skip the
download. Running an existing `GLOW/glow.exe` is not needed.

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
seasons and the two local solar times `12 h ∓ chi/15 h` (the equinox,
equator relation; dawn and dusk branches). With KCHEM=4 GLOW takes the electron
density above 200 km from this input, so the background must match the SZA:
an older background averaged over 00/06/12/18 LST and reused for every SZA put
a partly daytime ionosphere into the nightside columns. The night (chi=180°)
columns are 630.0 nm 94.9 R, 636.4 nm 30.0 R, 557.7 nm 12.8 R and 520 nm
2.62 R; at chi=0° 557.7 nm is 6.30 kR and 630 nm 2.53 kR. **These are results
for mean inputs, not a global climatological mean.**

The generator compiles patched sources in `airglow_data/build_glow/` and leaves
the originals unchanged. The patches (solzen removal, altitude grid, Earth
radius) are checked, and the build stops if a patch does not apply. The `data`
entry in this build is a relative symlink (`../../GLOW/data`), so keep `GLOW/`
and `airglow_data/` together in the project root or rebuild. The default NPZ/JSON
are overwritten on rerun.

## 6. Airglow step 2: combining PALACE nightglow components

### 6.1 Original archives and directories

The original distribution of the supplied PALACE v1.0 archives is
[Zenodo DOI 10.5281/zenodo.14064022](https://doi.org/10.5281/zenodo.14064022).
`PALACE/` already contains the three archives below and the extracted code and
data. **If an archive is missing**, download the missing `PALACE.zip`,
`test.zip` or `PMD.zip` from that Zenodo record into `PALACE/`, choosing
PALACE v1.0. The distribution can also be found through the
[code and data availability section of the PALACE paper](https://gmd.copernicus.org/articles/18/4353/2025/gmd-18-4353-2025-assets.html).
Existing files need not be downloaded or overwritten.

```text
PALACE/PALACE.zip
PALACE/test.zip
PALACE/PMD.zip
```

The model computation uses the code and data of `PALACE.zip`; `test.zip` holds
the distributed example for comparison. `PMD.zip` is auxiliary, but
`verify_palace.py` checks the CRC of all three archives, so **all three are
needed to run that verification.**

If they are not yet unpacked, run the code below with `python3`. Skip it if
`model/` and `reference_test/` already exist. To run it from a terminal, place
it between `python3 - <<'PY'` and a final `PY`.

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
are not implemented. See the local `AIRGLOW_COMPONENTS.md` for the latest
component assessment.

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
is the lowest supported altitude of the table (metadata
`support_min_altitude_km`, 80 km for the combined table). Below that altitude
results are flagged in `MODEL_VALID`; above the table top (1500 km) emission is
taken as zero without a flag. Add `--overwrite` to recompute existing files.

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

Use the current common interface. To avoid confusion with older CLI options of
the reference-data scripts, the final observation map is selected with
`--model` and the airglow component model with `--airglow-model`.

```bash
python3 sdrift_sky.py --model total --airglow-model combined --band VIS \
  --date 2026-10-19T00:00:00Z --orbit dawn-dusk --height 650 \
  --orbital-phase-deg 0 --min-elv 30 --sun-separation 90 --moon-separation 60 \
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

### 9.1 Sensitivity maps from Zodiacal + DGL + Airglow

After Section 8, the same `docs/data/` supports `model="sensitivity"`.
No separate sensitivity dataset, HST input file or additional Python package is
required. Keep `sensitivity.py` next to `sdrift_sky.py`, and include
`docs/assets/sensitivity.mjs` with the web assets.

The output is a **limiting surface brightness in AB mag/arcsec²**; a larger
limiting magnitude means a deeper detection limit. It is not a point-source
limiting magnitude. The model sums Zodiacal light, DGL and GLOW + PALACE
Airglow in linear intensity, then applies the Gaussian CCD equation with
source/sky shot noise, dark current, read noise and a sky-only systematic term.
It applies no HST correction or Earthshine contribution.

| Setting | Python / CLI | Default |
|---|---|---|
| Sky systematic fraction | `f_sys` / `--f-sys` | **0.001 = 0.1%** |
| Usable accumulated exposure | `exposure_hours` / `--exposure-hours` | **10 h** |
| Binning side length | `bin_side` / `--bin-side` | **5**, i.e. 5×5 detector pixels |
| Detection S/N | `snr` / `--snr` | 3 |
| Readout count | `n_reads` / `--n-reads` | 1 |

The instrument defaults from the reference notebook are diameter 30 cm,
detector pixel size 1.96″, read noise 1 e⁻/pixel/read, dark current
0.05 e⁻/pixel/s, and constant effective efficiency 50% for UV1/UV2/UV3 and
90% for VIS, with zero out-of-band leakage. The UV3 efficiency extends the UV
assumption of the notebook. A 5×5 bin covers **9.8″×9.8″**. Binning is independent
of HEALPix NSIDE. For multiple accumulated exposures, set `n_reads` to the actual
readout count; the default of one follows the notebook assumption.

Exposure means **usable accumulated time at the specified fixed sky/orbit
geometry**. It does not advance the orbit over those 10 hours or integrate
pixel-dependent visibility and background histories. The existing epoch/phase
option sets the satellite position at the selected UTC only.

Sensitivity requires Sun and Moon avoidance and `airglow_model="combined"`.
The current defaults are Sun separation ≥90°, Moon separation ≥60°, minimum
ELV 30°, H=650 km, dawn-dusk and north-referenced phase 0°. The separation
angles, ELV and orbit parameters remain adjustable. Directions rejected by
avoidance, ELV or Earth occultation are NaN. As in the other maps, sightlines
outside airglow model support **retain their partial calculation and are
hatched**, rather than being removed. Check `model_valid` or `valid` when
selecting samples for analysis; missing emission can make the partial
sensitivity optimistic.

```bash
mkdir -p output/fits/maps output/pdf/maps
python3 sdrift_sky.py --model sensitivity --band VIS \
  --date 2026-10-19T00:00:00Z --orbit dawn-dusk --height 650 \
  --orbital-phase-deg 0 --min-elv 30 --sun-separation 90 --moon-separation 60 \
  --f-sys 0.001 --exposure-hours 10 --bin-side 5 --snr 3 --n-reads 1 \
  --nside 128 --contour-step 0.2 \
  --fits output/fits/maps/sensitivity_VIS_nside128_20261019.fits.gz \
  --pdf output/pdf/maps/sensitivity_VIS_nside128_20261019.pdf
```

Use any of `UV1`, `UV2`, `UV3`, `VIS`, and add `--coords galactic` for Galactic
coordinates. The Python API uses the same settings:

```python
from sdrift_sky import SkyModel

sky = SkyModel()
r = sky.sky_map(
    "2026-10-19T00:00:00Z", nside=128, model="sensitivity", band="VIS",
    f_sys=0.001, exposure_hours=10, bin_side=5, snr=3, n_reads=1,
    orbit="dawn-dusk", height=650, orbital_phase_deg=0, min_elv_deg=30,
    sun_separation=90, moon_separation=60,
)
print(r["limiting_ab_mag_arcsec2"][r["valid"]])
print(r["background_mjy_sr"][r["valid"]])
```

In this mode, `mag` and `flux` are the **detection limit**, while
`background_mjy_sr` and the `zodiacal`, `dgl`, `airglow` arrays describe the
background. Sensitivity FITS differs from the 14-column background format:
the first two columns are `LIMIT_AB_MAG_ARCSEC2` and `LIMIT_MJY_SR`; the file
also contains `BACKGROUND_MJY_SR`, `SKY_COUNTS_E`, `SOURCE_COUNTS_E`, component
intensities and geometry flags. Headers record `FSYS`, `EXPHOURS`, `BINSIDE`,
`SNR`, `NREADS` and the instrument settings. Read by column name, not number.

For the web server above, select **Sensitivity · 한계 표면밝기** in the model
selector. The web `f_sys (%)` field uses **0.1**, equivalent to Python
`f_sys=0.001`. Change exposure, binning, S/N or readout count to recompute;
FITS.gz, PDF/print, PNG, Sun/Moon positions, separation boundaries and contours
are available as for the other maps. This option is in the web interface and
Python API/CLI, not the older desktop GUI.

For the full count equation, systematic floor, configurable Python instrument
parameters, notebook provenance and FITS details, see
[README_sensitivity.md](README_sensitivity.md). Changing sensitivity settings
requires only recalculating the final map; it does not require regenerating
the background dataset. Existing background model limitations remain.

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

Run the following tests included in the full handover (see `tests/README.md`).
Tests that compare with the JavaScript implementation call Node.js from Python,
so Node.js is required. The tests need the generated `docs/data/`, the
reference FITS and the airglow emission-rate NPZ files.

```bash
python3 -m unittest tests.test_dgl tests.test_airglow tests.test_airglow_components
python3 -m tests.test_web_sky
python3 -m unittest tests.test_airglow_orbit tests.test_sky_coordinates tests.test_python_interface tests.test_zodiacal_harmonics tests.test_orbit_gui_state
python3 -m unittest tests.test_sensitivity
python3 -m tests.test_moon_avoidance
node tests/test_web_dgl.mjs
node tests/test_web_moon.mjs
node tests/test_web_dates.mjs
node tests/test_web_worker_requests.mjs
node tests/test_web_app_cancel.mjs
```

The Python/JavaScript parity tests also cover the sun-synchronous orbit
elements, the Moon parallax and the zodiacal harmonics. Review PDFs and FITS
written by the tests go to `tmp/tests/`. The sensitivity tests additionally
check the CCD S/N equation, AB conversion, long-exposure systematic floor,
partial-airglow mask policy, Python/browser agreement across all four bands
and both coordinate systems, and the sensitivity FITS output. The Python test
runs `tests/test_sensitivity.mjs` using Node.js and temporary output files.

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
| Sensitivity f_sys, exposure, binning, S/N, readouts or instrument settings | no source regeneration; recompute the sensitivity map |

If `docs/data/`, which the common API reads, is not updated, rebuilding the NPZ
or reference FITS does not change the displayed values. Update the manifest and
its arrays together.

## 12. Archiving the handover and reproducibility

The provider keeps **the original inputs, generation code, validation code and
licenses** together, and hands over the full folder after it passes the input
check of Section 2.3 and the validation of Section 10. `.venv`, `__pycache__`,
Fortran executables and Cython binaries built on another computer are not used
as runtime dependencies in the new environment; compile the needed executables
again.

Run everything in a working copy that preserves the originals and the received
results. Section 2.4 moves the existing outputs to a separate folder so that the
rebuild does not reuse earlier results. Keep the new results together with the
code, inputs, commands, package versions and logs. Different package or compiler
versions can change the last digits or the gzip hash, so distinguish
bit-identical agreement from reproduction within numerical tolerance.

Preserve the GLOW academic research license, the GPLv3 of the PALACE code and
the CC BY 4.0 of its model data, and the attributions of the Zodiacal and DGL
sources. This procedure is completed locally; it performs no GitHub
commit/push and no update of the deployment folder.

## 13. Command block to run after preparation

From the root of the working copy after completing Sections 2.2–2.4, paste the
whole block below into a terminal. It does not archive old results or activate
the virtual environment. It performs **the actual recomputation**, so the
Zodiacal templates (about 2.5 h with 12 workers on the development machine,
longer with 4) and the LOS LUTs take a long time. Run times depend on the
computer.

By default it builds the reference data for the four bands, the two airglow
models, the H/SZA maps and PDFs, all shared data, the total maps for both orbits
and the validation. The LOS LUTs for the older viewer are built by default too.
If only the current web and common Python interface are needed, set
`SDRIFT_BUILD_LOS_LUT=0`. Adjust the number of local processes with
`SDRIFT_WORKERS`.

```bash
bash <<'BASH'
set -euo pipefail
SDRIFT_WORKERS=4
SDRIFT_BUILD_LOS_LUT=1
mkdir -p rebuild_logs output/pdf/airglow output/pdf/maps output/fits/maps
export MPLBACKEND=Agg

run_step() {
    local step_name="$1"
    shift
    "$@" 2>&1 | tee "rebuild_logs/${step_name}.log"
}

# Restore the pure-Python PALACE modules from the received archives.
python3 - <<'PY'
from pathlib import Path
from zipfile import ZipFile
folder = Path('PALACE')
for archive, source, destination in [
    ('PALACE.zip', 'PALACE', 'model'),
    ('test.zip', 'test', 'reference_test'),
]:
    target = folder / destination
    if not target.exists():
        staging = folder / 'unpacked'
        staging.mkdir(exist_ok=True)
        with ZipFile(folder / archive) as z:
            z.extractall(staging)
        (staging / source).rename(target)
PY

run_step 01_zodiacal python3 zodiacal_healpix.py --make-template \
    --nside 128 --band all --reference-date 2000-03-20T07:35:00 --epochs 12 \
    --workers "$SDRIFT_WORKERS" --no-plot
run_step 02_dgl python3 dgl_healpix.py --nside 128 --step-nm 0.5 --overwrite
run_step 03_glow python3 prepare_airglow_glow.py --f107 100 --ap 4 \
    --chi-step 1 --top-km 1500 --radius-km 6371
run_step 04_palace_verify python3 verify_palace.py
run_step 05_combined python3 prepare_airglow_palace.py --overwrite

for source in glow_selected_lines_mean glow_palace_night_mean; do
    run_step "06_${source}_maps" python3 airglow_healpix.py \
        --table "airglow_data/${source}.npz" --H 550 --sza 90 \
        --band UV1 UV2 UV3 VIS --nside 128 --workers "$SDRIFT_WORKERS" \
        --engine numpy --overwrite --pdf-dir output/pdf/airglow --contour-step 0.2
    if [ "$SDRIFT_BUILD_LOS_LUT" = 1 ]; then
        run_step "07_${source}_lut" python3 prepare_airglow_lut.py \
            --table "airglow_data/${source}.npz" --workers "$SDRIFT_WORKERS" --overwrite
        run_step "08_${source}_lut_check" python3 validate_airglow_lut.py \
            --table "airglow_data/${source}.npz" --nside 16
    fi
done
run_step 08_top2000_table python3 prepare_airglow_glow.py --top-km 2000 \
    --output airglow_data/glow_selected_lines_mean_top2000.npz
run_step 08_numerics python3 validate_airglow_numerics.py

run_step 09_export python3 export_web_data.py
run_step 10_physics_tests python3 -m unittest tests.test_dgl tests.test_airglow tests.test_airglow_components
run_step 11_web_test python3 -m tests.test_web_sky
run_step 12_interface_tests python3 -m unittest \
    tests.test_airglow_orbit tests.test_sky_coordinates tests.test_python_interface tests.test_zodiacal_harmonics tests.test_orbit_gui_state
run_step 13_sensitivity_test python3 -m unittest tests.test_sensitivity
run_step 13_moon_test python3 -m tests.test_moon_avoidance
run_step 13_web_dgl_test node tests/test_web_dgl.mjs
run_step 13_web_moon_test node tests/test_web_moon.mjs
run_step 13_web_dates_test node tests/test_web_dates.mjs
run_step 13_web_worker_test node tests/test_web_worker_requests.mjs
run_step 13_web_cancel_test node tests/test_web_app_cancel.mjs
if [ "$SDRIFT_BUILD_LOS_LUT" = 1 ]; then
    run_step 13_lut_tests python3 -m unittest tests.test_airglow_lut
fi

for orbit in dawn-dusk noon-midnight; do
    for band in UV1 UV2 UV3 VIS; do
        stem="total_${band}_nside128_20261019_${orbit}"
        run_step "14_${stem}" python3 sdrift_sky.py \
            --model total --airglow-model combined --date 2026-10-19T00:00:00Z \
            --orbit "$orbit" --height 650 --orbital-phase-deg 0 --band "$band" \
            --min-elv 30 --sun-separation 90 --moon-separation 60 \
            --nside 128 --contour-step 0.2 --engine numpy \
            --fits "output/fits/maps/${stem}.fits.gz" --pdf "output/pdf/maps/${stem}.pdf"
        stem="sensitivity_${band}_nside128_20261019_${orbit}"
        run_step "15_${stem}" python3 sdrift_sky.py \
            --model sensitivity --airglow-model combined --date 2026-10-19T00:00:00Z \
            --orbit "$orbit" --height 650 --orbital-phase-deg 0 --band "$band" \
            --min-elv 30 --sun-separation 90 --moon-separation 60 \
            --f-sys 0.001 --exposure-hours 10 --bin-side 5 --snr 3 --n-reads 1 \
            --nside 128 --contour-step 0.2 --engine numpy \
            --fits "output/fits/maps/${stem}.fits.gz" --pdf "output/pdf/maps/${stem}.pdf"
    done
done
printf '%s\n' 'Rebuild and tests finished. Also run the FITS/manifest/hash check of Section 10.1.'
BASH
```

The map computations in this block use NumPy explicitly to make cross-platform
reproduction easier. For speed, change `--engine numpy` to `--engine auto`.
Validation tolerances follow each test. Combining Numba parallelism with many
worker processes can oversubscribe CPU threads, so adjust to your resources.

The two orbit examples above are maps at a fixed orbital phase of 0° (the
northernmost point, `NORTH_V2`). They are not time averages over the orbit. For
Galactic output add `--coords galactic`. To check other components, choose
`--model zodiacal`, `dgl`, `zodiacal_dgl` or `airglow`. The block also makes
sensitivity maps for all four bands and both orbits with the defaults in
Section 9.1; these use the same fixed geometry, not an orbit-time average.

After validation, run the web server and GUI in separate terminals.

```bash
# Web: serves docs/ from the current root. Occupies the terminal; run separately.
python3 -m http.server 8765 --directory docs
```

```bash
# GUI: unset MPLBACKEND=Agg. Needs Tkinter and a graphical desktop.
env -u MPLBACKEND python3 sdrift_sky.py --gui
```

## 14. Output layout

Read this guide directly in the current folder: `README_first.md`.
`web_publish/` is not used in this handover and rebuild procedure.

| Path | Contents |
|---|---|
| `zodi_data/`, `dgl_data/` | original inputs and reference FITS; keep the paths the code reads |
| `airglow_data/` | emission-rate NPZ, optional LOS LUTs, reference airglow FITS and metadata |
| `docs/data/` | newly generated web and common Python runtime data |
| `output/pdf/zodiacal/` | Zodiacal PDFs produced during development |
| `output/pdf/airglow/` | airglow reference map PDFs |
| `output/pdf/maps/` | final background and sensitivity map PDFs for a given date and orbit |
| geometry PDFs in `output/pdf/` | generated explanatory figures and documents (see below) |
| `output/tex/` | TeX sources of the dawn-dusk geometry documents |
| `output/fits/maps/` | final background and sensitivity map FITS.gz for a given date and orbit |
| `rebuild_logs/` | logs and environment information of this rebuild |
| `rebuild_archive/` | copy of the received results kept before the rebuild |
| `tmp/` | screenshots, PNGs for PDF review, TeX build intermediates and test outputs |
| `tests/` | Python/JavaScript validation code; separate from the user programs |
| `figs/airglow_geometry_draft.pdf` | supplied input geometry sketch; distinct from generated output |

The geometry PDFs are regenerated as follows:

```bash
# output/pdf/airglow_geometry.pdf (needs reportlab)
python3 figs/draw_airglow_geometry.py
# output/pdf/dawn_dusk_geometry.pdf and output/pdf/dawn_dusk_airglow_geometry.pdf
# from output/tex/ (needs xelatex; intermediates in tmp/pdfs/dawn_dusk_build/)
bash figs/build_dawn_dusk_document.sh
```

See `output/README.md` for the locations of the organized results.
Some generation scripts write PDFs to the current directory by default
(`zodiacal_healpix.py` date maps without `--pdf`, for example), so specify the
result folder with `--pdf-dir` or `--pdf` as in the block above. The Zodiacal
template build skips PDFs with `--no-plot`, which does not affect the FITS output.

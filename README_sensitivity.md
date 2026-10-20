# S-DRIFT sensitivity map

`model="sensitivity"` computes the **limiting surface brightness
(AB mag/arcsec²)** that can be detected against the combined background of
zodiacal light + DGL + airglow. A larger limiting magnitude means that fainter
surface brightness can be detected. It is not a point-source limiting
magnitude. Python and the web page use the same Gaussian CCD equation.

## Defaults and the meaning of the exposure time

| Item | Python parameter | Default |
|---|---|---:|
| Fractional background systematic | `f_sys` | **0.001 = 0.1%** |
| Effective accumulated exposure | `exposure_hours` | **10 h** |
| Side of the spatial bin | `bin_side` | **5 (5×5 pixels)** |
| Detection S/N | `snr` | 3 |
| Number of readouts | `n_reads` | 1 |
| Aperture | `diameter_cm` | 30 cm |
| Detector pixel scale | `pixel_scale_arcsec` | 1.96″ |
| Read noise | `read_noise_e` | 1 e⁻ per pixel per read |
| Dark current | `dark_current_e_s` | 0.05 e⁻ per pixel per second |
| Constant in-band efficiency | `throughput` | UV1/UV2/UV3: 0.5, VIS: 0.9 |

The exposure time is the **effective accumulated exposure** given by the user.
The background at the specified UTC and orbit position is assumed to persist
for that time; the orbit is not advanced through the 10 hours, and the visible
time of each pixel is not integrated. The existing phase-evolution option only
places the satellite between the reference UTC and the requested UTC; it is not
an exposure integration.

The default 5×5 binning covers 25 detector pixels, i.e. **9.8″ × 9.8″**.
HEALPix NSIDE is only the sampling of the sky map and is unrelated to the
detector pixels or the binning. **Binning is assumed to be done in software
after readout:** each of the N = bin_side² pixels is read separately and
contributes RN² per read, so the read-noise variance is N · n_reads · RN².
On-chip binning (one read of the summed charge) would contribute RN² once per
read and give a slightly deeper limit. With the defaults the sky counts are
10⁴–10⁶ e⁻ in the bin, so the difference is below 0.001 mag for exposures of
1 h or more (about 0.002 mag for UV1 at 0.005 MJy/sr with 36 reads); it
matters only for very dark backgrounds, very short exposures or hundreds of
readouts.

The default single readout follows the reference notebook. For a stack of
several exposures, set `n_reads` to the actual number of readouts. The total
exposure time and the number of readouts are independent inputs.

The 50% efficiency for UV3 extends the notebook's UV efficiency in the same
way; it is an assumption. Out-of-band leakage is 0%, and no separate optical
loss, central obstruction or PSF aperture correction is applied. Interpret
`throughput` as the effective detection efficiency including QE.

## Web usage

Start a server in the project root:

```bash
python3 -m http.server 8000
```

Open `http://localhost:8000/docs/` and choose the **Sensitivity** entry
as the model (the web interface is in Korean). Besides the date, band,
coordinates, orbit, height, phase and ELV, you can change `f_sys (%)`, the
exposure time, the bin side, the S/N and the number of readouts. The web input
**0.1 means 0.1%**, the same as the Python input **0.001**. The aperture,
pixel scale, efficiency, read noise and dark current use the defaults in the
table above on the web page.

In Sensitivity mode, Sun and Moon avoidance are on and the airglow component is
fixed to GLOW + PALACE. The project's default observing conditions are Sun
separation ≥ 90°, Moon separation ≥ 60°, ELV ≥ 30°, height 650 km, dawn-dusk
orbit and phase 0° from the northernmost point; the angles and orbit settings
can be changed. No limit is shown for directions occulted by Earth or excluded
by the avoidance or ELV conditions.

**Outside the airglow model support, the partial value is shown and hatched,
as in the existing total-sky map.** The limit there can be optimistic because
part of the airglow is missing. FITS files store the model validity separately.
Sun/Moon positions and avoidance boundaries, the ELV boundary, the contour step
(default 0.2 mag), and FITS.gz and PNG/PDF export work as for the other maps.

## Python API

```python
from sdrift_sky import SkyModel

sky = SkyModel()
settings = dict(
    model="sensitivity", band="VIS",
    f_sys=0.001, exposure_hours=10, bin_side=5, snr=3, n_reads=1,
    height=650, orbit="noon-midnight", orbital_phase_deg=90,
    min_elv_deg=30, sun_separation=90, moon_separation=60,
)
date = "2026-10-17T00:00:00Z"

# The inputs 30, 10 are always RA, DEC in degrees.
p = sky.point(30, 10, date, **settings)
print(p["limiting_ab_mag_arcsec2"])
print(p["background_mjy_sr"], p["visible"], p["model_valid"])

# Equatorial RING HEALPix map; coords="galactic" is also supported.
r = sky.sky_map(date, nside=128, coords="equatorial", **settings)
sky.save_fits(r, "sensitivity_VIS_nside128_20261017.fits.gz")
fig = sky.plot(r, contour_step=0.2)
fig.savefig("sensitivity_VIS_nside128_20261017.png", dpi=160)
```

`evaluate(directions, date, model="sensitivity", ...)` is also supported;
`directions` are equatorial unit vectors of shape `(3,)` or `(N, 3)`, as in the
rest of the API. `sun_avoidance=False`, `moon_avoidance=False` and
`airglow_model="glow"` raise an error in this mode, because the map uses all
components and both avoidance conditions. Matplotlib figures and the CLI are
supported; the model list of the desktop GUI was not extended.

In sensitivity mode the common result fields `mag`, `ab_mag_arcsec2`, `flux`
and `intensity_mjy_sr` hold the **detection limit**. The background is in
`background_mjy_sr`. The aliases `limiting_ab_mag_arcsec2` and
`limiting_mjy_sr` are also provided. The components `zodiacal`, `dgl` and
`airglow` remain background intensities; their sum is the background, not the
limit.

- `visible`: the direction passes the Sun/Moon, ELV and Earth-occultation
  conditions.
- `model_valid`: the airglow calculation is inside the model support.
- `valid`: `visible & model_valid`.
- `mag` is NaN where `visible=False`. Where `visible=True` and
  `model_valid=False` the finite partial value remains; select `r['valid']`
  for statistics over complete-model samples only.

## When the background is already known

```python
from sensitivity import calculate_sensitivity

r = calculate_sensitivity(
    background_mjy_sr=0.2, band="VIS",
    f_sys=0.001, exposure_hours=10, bin_side=5,
)
print(r["limiting_ab_mag_arcsec2"])
print(r["sky_counts_e"], r["source_counts_e"])
```

This low-level function accepts a scalar or a NumPy array and keeps NaN
backgrounds as NaN. It takes no orbit or mask; those are handled by `SkyModel`
or by the caller. The instrument parameters can be changed with the names in
the table above.

## Command line

```bash
python3 sdrift_sky.py \
  --model sensitivity --band VIS --date 2026-10-17T00:00:00Z \
  --orbit noon-midnight --phase 90 --height 650 --min-elv 30 \
  --sun-separation 90 --moon-separation 60 \
  --f-sys 0.001 --exposure-hours 10 --bin-side 5 --snr 3 --n-reads 1 \
  --nside 128 --fits sensitivity_VIS.fits.gz --png sensitivity_VIS.png

python3 sdrift_sky.py \
  --model sensitivity --band VIS --date 2026-10-17T00:00:00Z \
  --orbit noon-midnight --phase 90 --ra 30 --dec 10 \
  --json sensitivity_point.json
```

`airglow_orbit.py` uses the same CLI. `--pdf`, `--show` and
`--coords galactic` are available. The CLI exposes the first five sensitivity
settings of the table; the other instrument parameters can be changed through
the Python API.

## Equations

The components are first added in linear intensity:

\[
I_\nu^{sky}=I_\nu^{zodi}+I_\nu^{DGL}+I_\nu^{airglow}.
\]

All three components are photon-counting band means in MJy/sr
(\(\langle I_\nu\rangle=\int I_\nu\,d\lambda/\lambda\,/\ln(\lambda_2/\lambda_1)\)),
so for a constant in-band efficiency \(q\) the sky electron rate of one detector
pixel follows exactly:

\[
r_{sky}=I_\nu^{sky}\frac{10^{-17}}{h}
\ln\frac{\lambda_2}{\lambda_1}\,A\,\Omega_{pix}\,q,
\quad A=\pi(D/2)^2,\quad
\Omega_{pix}=\left(\frac{p}{206264.80624709636}\right)^2.
\]

Here \(D\) is in cm, \(p\) in arcsec and \(h=6.62607015\times10^{-27}\) erg s.
This does not assume a flat spectrum; it uses the definition of the
**photon-counting band mean**. It must not be applied directly to other
weightings or to an efficiency that varies within the band.

With \(N=\mathrm{bin\_side}^2\) and \(t=3600\,\mathrm{exposure\_hours}\),

\[
B_{sky}=N r_{sky}t,\qquad
B=N(r_{sky}t+d t+n_{read}RN^2),\qquad
\sigma_{sys}=f_{sys}B_{sky}.
\]

For source electrons \(S\) and detection threshold \(k=\mathrm{snr}\), the
equation and its positive solution are

\[
k=\frac{S}{\sqrt{S+B+\sigma_{sys}^2}},\qquad
S=\frac{k^2}{2}\left[1+\sqrt{1+\frac{4(B+\sigma_{sys}^2)}{k^2}}\right].
\]

Dividing by the detected electrons for unit surface brightness gives the
limiting intensity, which is converted to AB mag/arcsec². The sky shot noise is
included separately from the systematic. The systematic applies only to the
sky, not to dark current, read noise or the source. It is assumed to be fully
correlated over time and over the bin, so long exposures approach
\(I_{lim}\to k f_{sys}I_{sky}\). As in the notebook, the finite-sample error of
a separate sky-estimation region is not added.

## Differences from the reference notebooks

The references are two Jupyter notebooks:

- `~/SKY_Background/K-DRIFT_space/sensitivity_2026_0923/plot_sensitivity_2026_0923.ipynb`:
  instrument settings, Gaussian CCD equation, AB conversion.
- `plot_sensitivity_2026_0923_sys.ipynb` in the same folder: the sky-only
  systematic term and the Gaussian closed form.

This implementation uses the **Gaussian method**; the Gehrels option of the
`_sys` notebook is not implemented. The notebook's HST sky spectrum is replaced
by the sum of the three S-DRIFT components. Earthshine, an empirical correction
for the difference from HST, and out-of-band leakage are not added. The
notebook sums a 1 Å grid; here the band integral for a flat throughput is done
analytically. Recomputing the notebook formula on its 1 Å grid with the same
backgrounds gives limits within 0.0004 mag of this implementation (grid end
points and older physical constants).

The limitations of the background models carry over to the sensitivity:
selected airglow components, no UV absorption, the night-layer assumptions, the
DGL UV extension, and so on. The map is not an exposure-schedule prediction
that follows the satellite motion.

## FITS and validation

The first two FITS columns are `LIMIT_AB_MAG_ARCSEC2` and `LIMIT_MJY_SR`;
`BACKGROUND_MJY_SR`, `SKY_COUNTS_E` and `SOURCE_COUNTS_E` are stored
separately. The background components and the geometry/validity flags are kept.
The Python and web files keep their own geometry column layouts, so read
columns by name, not by number. The header records `FSYS` (a fraction, 0.001),
`EXPHOURS` (hours), `BINSIDE`, `SNR`, `NREADS`, `DIAMCM`, `PIXSCALE`, `RN_E`,
`DARK_E_S`, `THRUPUT`, `NOISEMOD=GAUSSIAN`, `EXPMODE=USABLE_FIXED_SKY` and
`LEAKAGE=0`.

```bash
python3 -m unittest tests.test_sensitivity tests.test_python_interface -v
node tests/test_web_app_cancel.mjs
node tests/test_web_worker_requests.mjs
```

The tests cover the S/N obtained by substituting back into the CCD equation,
an independent AB zero-point integral, the long-exposure systematic floor,
trends with exposure, binning and readouts, the avoidance and partial-airglow
masks, agreement of the real web worker and Python maps for 4 bands × 2
coordinate systems, and the Python and web FITS columns and settings.

# 날짜와 이상적 궤도를 사용하는 airglow 지도

TLE가 없는 초기 설계 단계에서 두 원형 궤도를 비교합니다. Python과 웹 모두
`docs/data/`의 동일한 평균 발광률 표와 태양·달 위치 자료를 사용하며, 원본 GLOW나
PALACE를 설치할 필요가 없습니다. 출력 좌표는 기존 zodiacal/DGL 지도와 같은
지구 중심 GCRS 축의 RA, DEC입니다. 고도 기본값은 650 km, NSIDE는 128입니다.

## 실행

```sh
python3 -m pip install -r requirements-python.txt
# 선택: 현재 Python/NumPy와 호환되는 Numba를 설치하면 병렬 적분 사용
python3 -m pip install numba
python3 airglow_orbit.py --gui
```

`python3 airglow_orbit_gui.py`도 같은 창을 엽니다. 원본 계산 코드를 포함한 개발 폴더에서는
`python3 airglow_gui.py`도 새 뷰어를 엽니다. 그 폴더의 이전 H/SZA 기준 뷰어는
`python3 airglow_gui.py --reference`로 실행합니다. GitHub 배포에는 새 뷰어만 포함합니다.
GUI는 현재 UTC와 Zodiacal + DGL + Airglow 합산 모드로 시작합니다. Orbit, UTC, H, 북쪽 기준 위상 φ, band, 최소 ELV를 선택하고
Calculate map을 누릅니다. contour는 0.1-1.0 mag, 기본 0.5 mag입니다.
Tkinter와 그래픽 데스크톱이 필요합니다. `python3 -m tkinter`로 확인하고,
Tkinter가 없다면 Python 배포판에 맞는 Tcl/Tk 패키지를 설치합니다.
macOS Homebrew Python에서는 해당 Python 버전의 `python-tk@3.x`를 사용합니다.
`tkinter`를 pip로 설치하지 않습니다.

웹은 `docs/index.html`에서 모델을 **Airglow**로 선택합니다.

```sh
python3 -m http.server 8765 --directory docs
```

http://localhost:8765/ 로 접속합니다. 기존 GitHub Pages에서도 같은 UI를 사용합니다.
Web Worker에서 시선을 적분하며 계산을 중지할 수 있습니다. Python은 Numba가 있으면
병렬 적분하고, 없으면 NumPy로 실행합니다. `--engine numpy`로 선택할 수도 있습니다.

## 천정각, 태양·달 회피

천정각 θ는 위성의 지구 중심 반대 방향(천정)에서 잽니다.
0°는 천정, 90°는 수평, 180°는 지구 중심 방향입니다.
천정각은 포인터와 FITS의 참고 정보로 표시하며 최소·최대 천정각 입력과 필터는
사용하지 않습니다. 모든 시선에서 최소 ELV, 태양·달 회피, 지구 차폐 조건으로
관측 영역을 정합니다. 지구에 가리지 않는 하향 시선은 위성보다 낮은 대기를
통과하는 경로도 적분합니다.

**기본 관측 조건은 Sun separation >= 90°, Moon separation >= 60°, ELV >= 20°입니다.**
ELV(Earth limb elevation)는 구형 지구 표면의 가장자리보다 시선이 얼마나
위에 있는지를 나타내는 각도입니다. 고도 H와 천정각 θ에 대해
`ELV = 180° - asin(R/(R+H)) - θ`이며 R=6371 km입니다.
고도 650 km에서는 약 `ELV = 114.85° - θ`입니다.
GUI/웹의 **최소 ELV**를 조절하며 CLI는 `--min-elv 10`, Python API는
`min_elv_deg=10`처럼 지정합니다. 기본값은 20.0°이고 -90~180° 범위에서
입력합니다. 태양·달 회피 조건과 동시에 적용합니다.
음수를 지정해도 지구에 가려진 방향의 밝기는 계산하지 않습니다.
ELV >= 0°가 모델의 모든 대기층을 피한다는 뜻은 아닙니다.
최소 ELV 경계는 녹색 점선이며 포인터에 ELV와 천정각을 표시합니다.

Python과 웹 모두 태양·달을 표시하고 태양 경계는 빨간 점선, 달 경계는 주황 점선으로
그립니다. 회피각은 변경하거나 해제할 수 있습니다. 천정 위치는 파란색입니다.
회피·ELV 조건 밖은 밝은 무채색, 지구 차폐는 회색, 모델 범위 밖/광전자
cutoff는 빗금입니다. 빗금 영역의 밝기는 부분적인 계산 결과입니다.

## 날짜만으로 정해지지 않는 궤도상 위치

태양 단위벡터 s, 북극축 k로부터
`a = normalize(k - (k·s)s)`, `b = s × a`를 정의합니다.
선택한 날짜의 태양 방향을 사용하므로 춘분점에 고정된 지도가 아닙니다.

두 궤도 모두 **북쪽 기준 궤도 위상 φ**를 사용합니다. φ=0°는 북극축을
궤도면에 투영한 북쪽, φ=180°는 남쪽입니다. 지리적 북극 상공을 뜻하지 않습니다.

- Dawn-dusk: `r = (R+H)(cos(φ) a + sin(φ) b)`. 궤도면은 태양에 수직이고
  위성 SZA는 항상 90°입니다. φ=90°는 `b = s × a` 방향입니다.
- Noon-midnight: `r = (R+H)(cos(φ) a - sin(φ) s)`. 태양 방향을 포함하는
  궤도면입니다. φ=90°는 자정, φ=270°는 정오입니다.

이는 φ가 증가할 때 북쪽 → 자정 → 남쪽 → 정오로 이동한다는 정의이며,
이전 코드의 진행 방향을 유지합니다. 이전 noon-midnight의 u와는
`φ = (u - 90°) mod 360°` 관계이고, dawn-dusk에서는 이전 u와 같습니다.
**위상 기본값은 두 궤도 모두 북쪽인 0°**입니다. Noon-midnight의 야간 중심은
φ=90°를 선택합니다. 정오인 φ=270°에서는 태양 회피와 ELV 하한에 따라 관측 영역이 제한됩니다.

기본적으로 선택한 위상은 날짜를 바꿔도 유지합니다. 위성의 위치를 시간에 따라
이동시키려면 웹/GUI의 **기준 UTC부터 위상 자동 진행**을 선택하거나 CLI에
`--epoch`를 추가합니다. 기준 UTC의 위상이 입력한 `--orbital-phase-deg`이고,
`φ(t) = φ0 + sqrt(mu/(R+H)^3) (t-t0)`로 진행합니다(수식의 각도는 radian).
Python API 변수는 `orbital_phase_deg`이며 단위는 degree입니다. 기존 `--phase`와
API의 `phase=`도 별칭으로 받을 수 있지만 **새 북쪽 기준 정의로 해석**합니다.
이전 noon-midnight 값을 재현하려면 90°를 빼십시오. 기존 파일은 수정하지 않습니다.
R=6371 km, mu=398600.4418 km³/s²이며 H=650 km의 주기는 약 97.58분입니다.
기준 UTC 기본 예시 2026-10-17 00:00:00은 임의의 설계 기준입니다.
자동 진행 모드의 날짜 변경은 태양 방향, 달 방향, 궤도 위상에 모두 반영됩니다.

이 모델은 매 날짜 태양 방향을 따르는 **이상적인 기하 비교 모델**입니다.
고정 궤도경사각, 실제 LTAN, J2 섭동, SGP4/TLE 전파를 적용하지 않습니다.
실제 태양동기궤도와 두 이상적 평면이 정확히 같다는 뜻은 아닙니다.
궤도와 위상 convention에 대한 참고:
[NASA 궤도 개요](https://science.nasa.gov/earth/earth-observatory/catalog-of-earth-satellite-orbits/),
[Astropy 실제 위성 좌표/TLE 처리](https://docs.astropy.org/en/stable/coordinates/satellites.html).

## CLI와 Python API

```sh
# Dawn-dusk, 선택 UTC의 북쪽 위상, 기본 ELV 및 태양·달 회피
python3 airglow_orbit.py --kind airglow --date 2026-10-17T00:00:00Z \
  --orbit dawn-dusk --height 650 --orbital-phase-deg 0 --band VIS --nside 128 \
  --fits airglow_VIS_nside128_20261017_dawn-dusk.fits.gz \
  --png airglow_VIS_nside128_20261017_dawn-dusk.png

# Noon-midnight 자정; PDF 또는 --show로 화면 표시 가능
python3 airglow_orbit.py --kind airglow --date 2026-10-17T00:00:00Z \
  --orbit noon-midnight --orbital-phase-deg 90 --height 650 --band UV3 \
  --pdf airglow_UV3_nside128_20261017_noon-midnight.pdf

# 기준 위상부터 시간 진행, 단일 RA/DEC 질의
python3 airglow_orbit.py --kind airglow --date 2026-10-17T00:10:00Z \
  --epoch 2026-10-17T00:00:00Z --orbital-phase-deg 90 --orbit noon-midnight \
  --ra 30 --dec 10 --band VIS
```

`--sun-separation 100 --moon-separation 70`으로 회피각을 바꿉니다.
`--no-sun-avoidance`, `--no-moon-avoidance`로 각각 해제합니다.
단일 좌표 질의도 관측 조건을 적용합니다. 제외 방향은 mag/flux가 NaN이며,
`view_allowed`, `earth_occulted`, `model_valid`, 천정각과 separation을 함께 반환합니다.

```python
from airglow_orbit import AirglowOrbitModel
sky = AirglowOrbitModel()
r = sky.sky_map('2026-10-17T00:00:00Z', kind='airglow', orbit='noon-midnight',
                orbital_phase_deg=90, height=650, band='VIS', nside=128, min_elv_deg=20)
fig = sky.plot(r, contour_step=0.5)
fig.savefig('map.pdf')
sky.save_fits(r, 'map.fits.gz')
```

출력 FITS.gz는 RING/적도 HEALPix이며 AB_MAG_ARCSEC2, PHOTON_RADIANCE,
SUNLIT_RADIANCE, EARTH_OCCULTED, MODEL_VALID, MIN_ALTITUDE, PHOTO_CUTOFF,
ZENITH_ANGLE, VIEW_ALLOWED, ELV(deg)를 저장합니다. 제외 영역의 밝기는 NaN입니다.
VIEW_ALLOWED는 최소 ELV·태양·달 회피를 판정하고, 지구 차폐와 모델 유효성은
별도 열을 봅니다. 제외 방향은 적분을 생략하므로 다른 진단 열은 유효하지 않습니다.
UTC, 궤도 종류, 위상/기준 UTC, 고도, 위성 SZA, 회피각, 최소 ELV(ELVMIN),
태양·달·천정의 RA/DEC, 발광률 표 SHA-256 등을 헤더에 기록합니다.
위상 기준은 `PHREF=NORTH`, 정의 버전은 `PHCONV=NORTH_V1`로 명시합니다.
파일명에서는 궤도 위상을 `_phi...`로 표시합니다.

## 발광 모델의 범위

UV1=200-400, UV2=250-400, UV3=300-400, VIS=400-700 nm이며 constant
photon-counting throughput입니다. 고도와 **각 적분점의 국소 SZA**에서 발광률을
내삽하고, 지구 그림자를 적분점마다 판정하여 2 km 간격의 midpoint 적분을 합니다.
위성 SZA가 90°여도 시선상의 대기 SZA가 모두 90°라는 뜻은 아닙니다.
광자 radiance에서 band 평균 f_nu를 계산하여 AB mag/arcsec²로 변환합니다.

날짜는 기하학적 위치에 적용되고, 발광률은 기존 F10.7=100, Ap=4 등의 평균 표를
유지합니다. 기존 표의 누락 성분, UV 흡수 미포함, dawn/dusk 화학 이력 미포함,
PALACE 야간층의 전 지구 적용 가정도 유지됩니다. 달 위치는 기존 웹과 같은
지구 중심 위치이며 위성 시차와 산란 달빛은 포함하지 않습니다. GLOW 사용 조건은
`docs/data/glowlicense.txt`를 참조하십시오.

## Zodiacal light + DGL + Airglow 합산

웹 모델 목록에서 **Zodiacal light + DGL + Airglow**를 선택하면 Airglow와
동일한 궤도·위상·고도·ELV 입력을 사용할 수 있습니다. Python GUI에서는
Map을 **Zodiacal + DGL + Airglow**로 선택하거나 다음과 같이 합산 모드로 시작합니다.

```sh
python3 airglow_orbit.py --gui --kind total
python3 airglow_orbit.py --kind total --orbit noon-midnight \
  --orbital-phase-deg 90 --height 650 --date 2026-10-17T00:00:00Z \
  --band VIS --min-elv 5 --fits total_VIS_nside128_20261017.fits.gz \
  --png total_VIS_nside128_20261017.png
```

API에서도 `sky.sky_map(..., kind='total')`로 선택합니다. 기본값은 `kind='total'`로 세 성분을 합산합니다.
Airglow 단독은 `kind='airglow'` 또는 CLI의 `--kind airglow`로 선택합니다. `--model combined/glow`는 합산 여부와 별개로
Airglow의 발광 성분 모델을 선택합니다.

세 성분은 같은 band의 **MJy/sr로 통일하여 합한 뒤** AB mag/arcsec²로 변환합니다.
광자 radiance I_gamma의 변환은
`I_airglow[MJy/sr] = 6.62607015e-10 * I_gamma / ln(lambda_max/lambda_min)`이며
constant photon-counting throughput 가정과 기존 Airglow AB 정의를 사용합니다.
등급 자체를 더하지 않습니다. Zodiacal은 날짜에 따라 회전하고 DGL은 RA/DEC에
고정되며, Airglow는 선택한 궤도의 시선을 따라 적분합니다.

세 성분 모두 동일한 Sun/Moon avoidance, 최소 ELV, 지구 차폐를 적용합니다.
Airglow 모델이 불완전한 방향은 합산 지도에도 빗금으로 표시합니다.
총합 FITS의 두 번째 열은 `INTENSITY_MJY_SR`이고 별도로
`ZODIACAL_MJY_SR`, `DGL_MJY_SR`, `AIRGLOW_MJY_SR`,
`AIRGLOW_PHOTON_RADIANCE`를 저장합니다. 모델 유효성 등 기존 진단 열도 유지합니다.
총합 결과의 API `flux` 단위는 MJy/sr이고 Airglow 단독의 `flux`는 photons/s/cm²/sr입니다.
`airglow_photon_radiance`에 합산 전 Airglow 광자 radiance도 보관합니다.

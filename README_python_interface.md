# S-DRIFT Python 인터페이스와 사용법

작성 기준: 2026-10-19, 공통 Python API v2.

Zodiacal light, DGL, Airglow와 합산 밝기를 **하나의 `SkyModel` 인터페이스**로
계산합니다. 모드가 달라도 함수 이름·설정 인자·intensity 단위·관측 조건을
동일하게 사용합니다. 이 문서는 Python API, 명령행(CLI), GUI와 저장 파일을 설명합니다.

## 1. 가장 간단한 사용법

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

`point(0, 60, 날짜, ...)`에서 **0은 RA=0°, 60은 DEC=+60°**입니다.
RA의 단위는 시간이 아닌 degree이며 RA=0°는 0시간에 해당합니다.
날짜 끝의 `Z`는 UTC입니다. 위 예제는 해당 UTC, 해당 하늘 방향에서의
VIS 밴드 Zodiacal light + DGL 밝기를 계산합니다.

기본으로 궤도·ELV·태양/달 회피 조건을 적용합니다. `visible=False`인 방향의
밝기와 intensity는 `NaN`입니다. 따라서 밝기와 함께 관측 가능 여부를 확인합니다.

모드를 바꿀 때는 **`model`만 바꾸면 됩니다.**

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

## 2. 필요한 파일과 설치

GitHub 저장소를 내려받은 뒤 다음 파일 구조를 유지합니다.

```text
SDRIFT_sky/
├── sdrift_sky.py               # 공통 SkyModel, CLI
├── airglow_orbit.py            # 수치 계산 엔진, 호환 CLI/클래스 이름
├── airglow_orbit_gui.py        # 공통 Python GUI
├── sky_coordinates.py         # Equatorial / Galactic 회전
├── requirements-python.txt
└── docs/data/
    ├── manifest.json          # 파일명, 단위, 배열 크기, 출처와 해시
    └── ...                    # 배포된 압축 데이터 전체
```

Python 파일 하나만 분리해서 복사하지 말고, 위 파일과 `docs/data/`를 함께 둡니다.
기존 대용량 원본 FITS/NPZ, GLOW/PALACE 소스나 컴파일러는 배포 데이터의 조회에
필요하지 않습니다. 실행 중 인터넷 연결도 필요하지 않습니다.

```sh
python3 -m pip install -r requirements-python.txt
# 선택: Airglow 병렬 적분을 사용하려면 호환되는 Numba 설치
python3 -m pip install numba
```

Python 3.10 이상을 대상으로 합니다. NumPy, healpy, Astropy, Matplotlib가 필요하고
Numba는 선택 사항입니다. 설치된 Python·NumPy와 호환되는 Numba를 사용합니다.
GUI에는 그래픽 데스크톱, Tkinter와 해당 Python에 맞는 Tcl/Tk가 필요합니다.

```sh
python3 -m tkinter
```

위 명령으로 Tkinter 창이 열리는지 확인합니다. Tkinter는 Python 배포판의
패키지로 설치하며 `pip install tkinter`를 사용하지 않습니다.
화면 없는 서버에서는 GUI나 `--show` 대신 파일 저장을 사용합니다.

다른 데이터 폴더를 사용하려면 다음과 같이 지정합니다.

```python
from sdrift_sky import SkyModel
sky = SkyModel(data_dir="/absolute/path/to/docs/data")
```

데이터는 처음 사용할 때 읽고 크기와 SHA-256을 검증한 뒤 인스턴스에 캐시합니다.
여러 방향이나 날짜를 계산할 때 `sky` 인스턴스를 재사용합니다.

## 3. 공통 모드와 기본값

### 3.1 `model`: 합산할 하늘 성분

| 값 | 계산 성분 |
|---|---|
| `"zodiacal"` | Zodiacal light |
| `"dgl"` | Diffuse Galactic Light |
| `"zodiacal_dgl"` | Zodiacal light + DGL |
| `"airglow"` | Airglow |
| `"total"` | Zodiacal light + DGL + Airglow; 기본값 |

모든 모드가 같은 궤도와 관측 조건을 사용합니다. Zodiacal/DGL 모드에서도
관측 영역을 제한하지만, 선택하지 않은 Airglow를 계산해 더하지는 않습니다.

### 3.2 `airglow_model`: Airglow 내부 성분 선택

| 값 | 의미 |
|---|---|
| `"combined"` | GLOW + 가정한 PALACE 야간 발광층; 기본값 |
| `"glow"` | GLOW의 선택된 발광선 |

`model`과 `airglow_model`은 별개의 설정입니다.

```python
from sdrift_sky import SkyModel
sky = SkyModel()
p = sky.point(
    30, 10, "2026-10-17T00:00:00Z",
    model="total", airglow_model="glow",
    orbit="noon-midnight", orbital_phase_deg=90,
)
```

`model="zodiacal"`, `"dgl"`, `"zodiacal_dgl"`에서는 `airglow_model`을 합산에
사용하지 않습니다. GUI에서도 이 경우 Airglow 성분 선택을 비활성화합니다.

### 3.3 공통 입력 표

`point()`, `evaluate()`, `sky_map()`에 다음 계산 설정을 전달할 수 있습니다.

| Python 인자 | 기본값 | 단위·허용값·의미 |
|---|---|---|
| `model` | `"total"` | 위 다섯 하늘 모드 |
| `airglow_model` | `"combined"` | `combined`, `glow` |
| `band` | `"VIS"` | `UV1`, `UV2`, `UV3`, `VIS` |
| `orbit` | `"dawn-dusk"` | `dawn-dusk`, `noon-midnight` |
| `height` | `650.0` | km; `0 < height < 1500` |
| `orbital_phase_deg` | `0.0` | 북쪽 기준 위상, degree; 360° 주기로 해석 |
| `epoch` | `None` | 위상의 기준 UTC; `None`이면 위상 고정 |
| `min_elv_deg` | `20.0` | degree; −90~180 |
| `sun_avoidance` | `True` | 태양 회피 적용 여부 |
| `sun_separation` | `90.0` | degree; 0~180 |
| `moon_avoidance` | `True` | 달 회피 적용 여부 |
| `moon_separation` | `60.0` | degree; 0~180 |
| `engine` | `"auto"` | `auto`, `numba`, `numpy` |
| `step_km` | `2.0` | Airglow 적분의 최대 구간 길이, km; 양수; API에서 지정 |

`coords`는 기본 `"equatorial"`이며 `"galactic"`도 선택할 수 있습니다.
좌표 입력 및 표시 규칙은 5절을 참조하십시오.
`nside`는 `sky_map()`만, `contour_step`은 `plot()`만의 인자입니다.
`data_dir`는 `SkyModel()`을 생성할 때 지정합니다.

### 3.4 밴드와 밝기 단위

| Band | 파장 범위 |
|---|---|
| `UV1` | 200–400 nm = 0.2–0.4 μm |
| `UV2` | 250–400 nm = 0.25–0.4 μm |
| `UV3` | 300–400 nm = 0.3–0.4 μm |
| `VIS` | 400–700 nm = 0.4–0.7 μm |

모든 밴드는 일정한 photon-counting throughput을 가정합니다.
**모든 모드의 intensity는 MJy/sr**, 밝기는 **AB mag/arcsec²**로 통일했습니다.
Airglow 단독에서도 같은 단위를 사용합니다.
별도의 `photon_radiance` 또는 `airglow_photon_radiance`는 Airglow의 광자 단위 값입니다.
선택하지 않은 성분은 관측 가능 방향에서 0입니다.

합산은 같은 밴드의 intensity를 더한 뒤 AB 등급으로 변환합니다.
AB 등급을 직접 더하지 않습니다. 큰 AB 등급일수록 더 어둡습니다.

## 4. 날짜, 궤도와 관측 조건

### 4.1 UTC

날짜는 ISO 형식의 문자열로 전달합니다.

```text
2026-10-19T00:00:00Z
2026-10-19T09:00:00+09:00
```

위 두 시각은 같습니다. 시간대가 없는 문자열도 받지만 UTC로 해석하므로
`Z` 또는 시간대 오프셋을 명시하는 편이 안전합니다.
저장된 태양·달 위치 데이터의 사용 범위는 1900–2100년입니다.
API의 `date`는 필수이고 CLI/GUI에서 생략하면 현재 UTC를 사용합니다.

### 4.2 북쪽 기준 궤도 위상 φ

지구 중심에서 북쪽 기준으로 잰 궤도상 위치입니다. 두 궤도 모두 0°는 북쪽,
180°는 남쪽이며, 북쪽은 지구 북극축을 궤도면에 투영한 방향입니다.

| 위상 | Dawn–dusk | Noon–midnight |
|---|---|---|
| 0° | 북쪽 | 북쪽 |
| 90° | 태양 방향 × 북쪽 방향 | 자정 |
| 180° | 남쪽 | 남쪽 |
| 270° | 90° 방향의 반대 | 정오 |

Dawn–dusk 궤도면은 태양에 수직이고 위성 SZA는 90°입니다.
Noon–midnight 궤도면은 태양 방향을 포함합니다.
현재 모델은 TLE 없는 설계 단계의 원형 궤도 기하이며 실제 J2/SGP4 전파가 아닙니다.

### 4.3 위상 고정과 시간 진행

`epoch=None`이면 날짜를 바꿔도 입력한 위상을 유지합니다.
태양·달 위치와 천구상의 배치는 날짜에 맞게 바뀝니다.
`epoch`를 주면 그 UTC에 입력 위상에 있었다고 가정하고 공전 주기에 따라 진행합니다.

```python
from sdrift_sky import SkyModel
sky = SkyModel()
p = sky.point(
    30, 10, "2026-10-17T00:10:00Z",
    model="total", orbit="noon-midnight", height=650,
    orbital_phase_deg=90, epoch="2026-10-17T00:00:00Z",
)
print(p["orbital_phase_deg"])  # 약 126.893°
```

고도 650 km에서 공전 주기는 약 97.58분입니다. 기준 UTC는 설계상 가정입니다.

### 4.4 ELV와 회피각

천정각 θ는 위성의 천정에서 시선까지의 각도입니다.
0°는 천정, 90°는 수평, 180°는 천저입니다. 최소·최대 천정각 필터는 없습니다.

구형 지구 반지름 R=6371 km에서:

```text
ELV = 180° − asin[R / (R + H)] − θ
```

ELV=0°는 지표면에 접하는 방향이고 양수는 그보다 바깥쪽 하늘입니다.
고도 650 km에서는 `ELV ≈ 114.85° − θ`이며 기본 ELV≥20°는 θ≤94.85°입니다.

기본 관측 조건은 다음 세 가지를 동시에 만족해야 합니다.

- ELV ≥ 20°
- 태양과의 separation ≥ 90°
- 달과의 separation ≥ 60°

여기에 지구 차폐를 제외합니다. 회피 체크를 해제해도 지구 차폐는 해제되지 않습니다.
ELV를 음수로 낮춰도 지구에 가려진 방향의 밝기는 반환하지 않습니다.

## 5. 좌표계와 방향 입력

### 5.1 `point(ra, dec, date, ...)`

**입력 RA, DEC는 항상 적도 좌표 degree**입니다.
`coords="galactic"`을 지정해도 `point()`의 첫 두 숫자를 l, b로 해석하지 않습니다.
RA는 360° 주기로 정규화하고 DEC 범위는 −90~90°입니다.

### 5.2 `sky_map(..., coords=...)`

- `equatorial`: RA, DEC HEALPix 지도
- `galactic`: 은하 경도 l, 위도 b HEALPix 지도

지도·태양·달·천정·회피 경계를 같은 좌표계로 표시합니다.
FITS도 선택한 좌표계의 픽셀 격자로 저장합니다.
`coords`는 `point()`와 `evaluate()`에서도 결과의 좌표계 메타데이터를 지정할 수 있지만,
이들 함수의 입력 RA/DEC 또는 방향 벡터의 정의를 바꾸지는 않습니다.

### 5.3 l, b 방향을 단일 좌표로 질의하는 방법

```python
from sdrift_sky import SkyModel, vector, coordinates
from sky_coordinates import to_equatorial

sky = SkyModel()
l, b = 120.0, 30.0
ra, dec = coordinates(to_equatorial(vector(l, b), "galactic"))
p = sky.point(ra, dec, "2026-10-17T00:00:00Z", model="total")
print(ra, dec, p["ab_mag_arcsec2"], p["visible"])
```

기존 지구 중심 태양·달 벡터의 축을 고정 회전 행렬로 변환합니다.
위성 시차나 추가 광행차를 계산하는 기능은 아닙니다.

## 6. 단일 좌표 반환값: `point()`

반환값은 Python `dict`이며, 각 방향의 값은 스칼라입니다.

| 키 | 단위·의미 |
|---|---|
| `ra_deg`, `dec_deg` | 입력 방향의 적도 좌표 degree |
| `utc`, `band`, `model`, `airglow_model`, `coords` | 계산 설정 |
| `ab_mag_arcsec2` | AB mag/arcsec²; 제외 방향은 NaN |
| `intensity_mjy_sr` | 선택 성분의 총 intensity, MJy/sr |
| `zodiacal_mjy_sr`, `dgl_mjy_sr`, `airglow_mjy_sr` | 각 성분, MJy/sr |
| `photon_radiance` | Airglow만의 photons/s/cm²/sr; 비선택 시 0 |
| `view_allowed` | ELV와 태양·달 회피 조건 충족 |
| `earth_occulted` | 지구 차폐 여부 |
| `visible` | `view_allowed`이며 지구에 가리지 않음 |
| `model_valid` | 현재 수치 모델의 지원 범위·cutoff 판정 |
| `valid` | `visible and model_valid` |
| `display_ab_mag_arcsec2` | 관측 가능하면 밝기, 아니면 `None` |
| `orbit`, `height_km`, `orbital_phase_deg` | 궤도 종류, km, 해당 UTC의 위상 degree |
| `sza_deg`, `zenith_deg`, `elv_deg`, `minimum_elv_deg` | 태양천정각·천정각·ELV·하한, degree |
| `sun_ra_dec_deg`, `moon_ra_dec_deg` | 태양·달의 지구 중심 적도 좌표 `(RA, DEC)` |
| `sun_separation_deg`, `moon_separation_deg` | 시선과 태양·달 사이 각도, degree |
| `sun_avoidance`, `moon_avoidance` | 회피 적용 여부 |
| `minimum_sun_separation_deg`, `minimum_moon_separation_deg` | 설정 회피각, degree |

`visible=True`라도 `model_valid=False`일 수 있습니다. 이 경우 관측 기하는 가능하지만
밝기는 부분적인 계산 결과입니다. `model_valid=True`도 누락 발광 성분까지 모두
계산했다는 뜻은 아닙니다. 제외된 방향의 모델 판정은 해석하지 마십시오.

## 7. 전천 지도: `sky_map()`

```python
from pathlib import Path
from sdrift_sky import SkyModel

sky = SkyModel()
result = sky.sky_map(
    "2026-10-17T00:00:00Z",
    model="total", airglow_model="combined", band="VIS",
    nside=128, coords="galactic",
    orbit="noon-midnight", height=650, orbital_phase_deg=90,
    min_elv_deg=20, sun_separation=90, moon_separation=60,
)
Path("output/python").mkdir(parents=True, exist_ok=True)
sky.save_fits(result, "output/python/total_VIS_nside128_20261017_galactic.fits.gz")
fig = sky.plot(result, contour_step=0.2)
fig.savefig("output/python/total_galactic.png", dpi=160)
```

계산 설정은 `point()`와 같고, `nside`와 지도 좌표계를 추가로 지정합니다.

### 7.1 해상도와 픽셀 순서

`nside`는 양의 2의 거듭제곱이며 기본 128입니다. GUI에서는 16/32/64/128을
선택합니다. API와 CLI는 다른 유효한 NSIDE도 받을 수 있습니다.

```text
픽셀 수 N = 12 × nside²
NSIDE 16  →   3,072개
NSIDE 32  →  12,288개
NSIDE 64  →  49,152개
NSIDE 128 → 196,608개
```

배열 순서는 **RING HEALPix**입니다. 높은 NSIDE는 계산량을 늘리며, 입력
Zodiacal/DGL 템플릿 자체의 NSIDE 128보다 더 세밀한 관측 정보를 만들어 주지는 않습니다.
은하 좌표 지도는 은하 좌표 픽셀 중심의 실제 시선을 평가합니다.

### 7.2 지도 딕셔너리의 주요 키

방향별 배열의 크기는 `(N,)`입니다. 다음 반환 구조는 모든 `model`에서 같습니다.

| 키 | 내용 |
|---|---|
| `flux`, `intensity_mjy_sr` | 같은 총 intensity 배열, MJy/sr |
| `flux_unit` | 항상 `"MJy/sr"` |
| `mag`, `ab_mag_arcsec2` | 같은 AB mag/arcsec² 배열 |
| `zodiacal`, `dgl`, `airglow` | 각 성분 배열, MJy/sr |
| `airglow_photon_radiance` | Airglow 광자 radiance 배열 |
| `lit` | Airglow의 햇빛 받는 경로 부분 광자 radiance |
| `view`, `view_allowed` | ELV·태양·달 회피 조건의 동일 마스크 |
| `visible` | 회피 조건을 만족하고 지구에 가리지 않는 마스크 |
| `earth`, `earth_occulted` | 지구 차폐의 동일 마스크 |
| `model_valid` | 수치 모델 지원 범위 판정 |
| `valid` | `visible & model_valid` |
| `theta`, `elv` | 천정각과 ELV 배열, degree |
| `hmin` | 전방 시선 최소 고도, km; 지구 차폐 시 0 |
| `cutoff` | Airglow 광전자 cutoff 구간 포함 여부 |
| `nside`, `coords`, `model`, `band`, `utc` | 지도 설정 |
| `airglow_model`, `orbit`, `height`, `epoch`, `min_elv_deg` | 계산 설정 |
| `state` | 위성 기하; 아래 표 참조 |
| `sun`, `moon` | 적도 축의 단위벡터, 각각 `(3,)` |
| `engine` | 실제 계산 경로: `numba`, `numpy`, `map interpolation` |

`state`의 주요 키는 다음과 같습니다.

| 키 | 의미 |
|---|---|
| `zenith` | 해당 UTC 위성의 천정 단위벡터, 적도 축 |
| `orbital_phase_deg` | 기준 UTC로부터 진행시킨 후의 위상, degree |
| `sza` | 위성 위치의 태양천정각, degree |
| `period` | 공전 주기, **초** |
| `limb_zenith_deg` | 지표면 접선 방향의 천정각, degree |
| `axes` | 수치 적분의 기준 축에서 실제 적도 축으로의 회전행렬 `(3,3)` |

반환 딕셔너리는 변형하지 않고 `plot()`과 `save_fits()`에 전달합니다.
`_engine_result` 등 밑줄로 시작하는 항목은 내부 구현용입니다.

### 7.3 관측 가능한 방향의 통계

```python
import numpy as np
from sdrift_sky import SkyModel

sky = SkyModel()
r = sky.sky_map("2026-10-17T00:00:00Z", model="total", nside=16)
use = r["visible"] & r["model_valid"] & np.isfinite(r["mag"])
if use.any():
    print("선택한 픽셀 수:", np.count_nonzero(use))
    print("AB 등급 범위:", r["mag"][use].min(), r["mag"][use].max())
else:
    print("조건을 만족하는 유효 픽셀이 없습니다.")
```

등급을 단순 평균한 값은 평균 intensity의 등급과 다릅니다.
밝기의 평균이 필요하면 `flux`를 평균한 뒤 `ab_magnitude()`로 변환합니다.

## 8. 여러 시선의 직접 계산: `evaluate()`

`directions`는 **적도 축의 단위벡터** `(3,)` 또는 `(N,3)`입니다.
반환하는 방향별 값은 단일 벡터 입력도 `(1,)` 배열입니다.
RA, DEC 숫자를 그대로 벡터 자리에 넣지 않습니다.

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

`evaluate()` 결과는 HEALPix 지도가 아니므로 바로 `plot()`이나 `save_fits()`에
전달하지 않습니다. 전천 결과 저장에는 `sky_map()`을 사용합니다.

## 9. 지도 표시와 저장

### 9.1 `plot(result, target=None, contour_step=0.2, figure=None)`

`sky_map()` 결과를 받아 Matplotlib `Figure`를 반환합니다.

| 인자 | 의미 |
|---|---|
| `result` | `sky_map()` 반환 딕셔너리 |
| `target` | 선택 사항; 표시할 적도 좌표 `(RA, DEC)`, degree |
| `contour_step` | 기본 0.2 mag; 0.1~1.0 mag 중 0.1 간격 |
| `figure` | 기존 Matplotlib Figure에 그릴 때 지정; 생략 시 새 Figure |

이미지를 저장하려면 `fig.savefig("map.png")` 또는 `fig.savefig("map.pdf")`,
화면을 띄우려면 `matplotlib.pyplot.show()`를 호출합니다.
화면 없는 환경은 Python을 실행하기 전에 `MPLBACKEND=Agg`를 지정하거나,
`matplotlib.pyplot`을 import하기 전에 `matplotlib.use("Agg")`를 호출합니다.

색상 범위는 유효한 밝기 범위에서 자동 설정하며, 큰 등급은 어두운 색입니다.
태양은 빨간 별과 빨간 점선, 달은 노란/주황 원과 주황 점선,
최소 ELV는 녹색 점선, 천정은 파란 점으로 표시합니다.
회색은 지구 차폐, 옅은 영역은 제외 방향, 빗금은 Airglow의 부분적인 계산 방향입니다.

### 9.2 `save_fits(result, filename)`

파일명은 반드시 `.fits.gz`로 끝나야 합니다. 상위 폴더가 없으면 만들고,
같은 이름의 기존 파일은 덮어씁니다. API v2는 모든 모드에서 다음 **14개 열**을
같은 순서와 단위로 저장합니다.

| 열 이름 | 단위 |
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

열은 float32로 저장하므로 0/1 열은 읽은 뒤 필요에 따라 bool로 변환합니다.
`COORDSYS=C`는 적도 좌표, `G`는 은하 좌표이며 `ORDERING=RING`입니다.
헤더에는 `APIVER=2`, `MODEL`, `AGMODEL`, `FLUXUNIT=MJy/sr`, UTC, 위상,
고도, ELV·회피 조건, 데이터 출처와 해시 등을 기록합니다.

```python
from astropy.io import fits

with fits.open("output/python/total_VIS_nside128_20261017_galactic.fits.gz") as hd:
    header = hd[1].header
    table = hd[1].data
    mag = table["AB_MAG_ARCSEC2"].ravel()
    visible = table["VISIBLE"].ravel().astype(bool)
    print(header["COORDSYS"], header["APIVER"], mag[visible].shape)
```

이것은 **공통 Python API v2의 스키마**입니다. 기존 파일이나 웹 다운로드에는
과거의 열 구성이 있을 수 있으므로 열 번호보다 열 이름을 사용합니다.
태양·달·천정의 RA/DEC 헤더 값은 지도 좌표계와 무관하게 적도 좌표입니다.

### 9.3 NaN, infinity와 마스크

- `NaN`: 회피·ELV 조건 또는 지구 차폐로 제외된 방향 등 값이 없는 경우
- `+inf` 등급: 계산된 intensity가 0인 경우
- `visible`: 관측 기하 판정
- `model_valid`: 해당 수치 모델의 지원 범위 판정
- `valid`: 두 판정을 모두 만족

NaN과 infinity를 JSON으로 저장하려면 `None`으로 바꿉니다.
CLI의 단일 좌표 `--json` 출력은 이 변환을 자동으로 수행합니다.

## 10. 명령행: 두 파일에서 같은 옵션 사용

다음 두 명령은 **동일한 CLI**를 실행합니다.

```sh
python3 sdrift_sky.py --help
python3 airglow_orbit.py --help
```

### 10.1 단일 좌표

```sh
python3 sdrift_sky.py --ra 0 --dec 60 \
  --date 2026-10-13T12:00:00Z --model zodiacal_dgl --band VIS
```

`--ra`, `--dec`는 함께 지정합니다. 단일 좌표만 요청하면 지도는 만들지 않습니다.
단일 좌표와 `--png`/`--fits` 등을 같이 지정하면 좌표 질의 후 지도도 만듭니다.
`--coords galactic`을 주더라도 `--ra`, `--dec`는 적도 좌표입니다.

### 10.2 은하 좌표 총합 지도

```sh
python3 sdrift_sky.py --model total --airglow-model combined \
  --date 2026-10-17T00:00:00Z --band VIS --coords galactic --nside 128 \
  --orbit noon-midnight --orbital-phase-deg 90 --height 650 --min-elv 20 \
  --contour-step 0.2 --fits output/python/total_galactic.fits.gz \
  --png output/python/total_galactic.png
```

### 10.3 같은 조건에서 성분만 변경

```sh
python3 sdrift_sky.py --model dgl --date 2026-10-17T00:00:00Z \
  --orbit noon-midnight --orbital-phase-deg 90 --png output/python/dgl.png
python3 sdrift_sky.py --model airglow --airglow-model glow \
  --date 2026-10-17T00:00:00Z --orbit noon-midnight --orbital-phase-deg 90 \
  --png output/python/airglow.png
```

### 10.4 주요 CLI 옵션과 API 대응

| CLI | Python API / 의미 |
|---|---|
| `--model` | `model` |
| `--airglow-model` | `airglow_model` |
| `--date`, `--utc` | `date`; CLI 기본 현재 UTC |
| `--ra`, `--dec` | `point()`의 적도 좌표 |
| `--orbit` | `orbit` |
| `--height`, `--H` | `height` |
| `--orbital-phase-deg`, `--phase` | `orbital_phase_deg` |
| `--epoch` | `epoch`; 생략 시 위상 고정 |
| `--min-elv`, `--elv-min` | `min_elv_deg` |
| `--sun-separation`, `--moon-separation` | 각 회피각 |
| `--no-sun-avoidance`, `--no-moon-avoidance` | 대응 bool을 False로 설정 |
| `--band`, `--coords`, `--nside` | 같은 이름의 지도 설정 |
| `--engine` | `engine` |
| `--data-dir` | 생성자의 `data_dir` |
| `--contour-step` | `plot()`의 `contour_step` |
| `--show` | 지도 창 열기 |
| `--png`, `--pdf`, `--fits` | 각 출력 경로 |
| `--json` | 단일 좌표 dict 저장; `--ra`, `--dec` 필요 |
| `--gui` | 공통 GUI 열기 |

출력 옵션과 좌표를 모두 생략하면 지도를 계산하고 상태를 출력하지만
파일 저장이나 창 표시는 하지 않습니다. 원하는 출력 옵션을 추가하십시오.
`step_km`은 현재 API 인자이며 별도 CLI 옵션은 없습니다.

## 11. GUI

```sh
python3 sdrift_sky.py --gui
python3 airglow_orbit.py --gui
```

두 명령은 같은 창을 엽니다. 다음과 같이 초기 설정을 전달할 수도 있습니다.

```sh
python3 sdrift_sky.py --gui --model dgl --coords galactic \
  --orbit noon-midnight --orbital-phase-deg 90 --height 650 --min-elv 25
```

| 화면 항목 | 의미 |
|---|---|
| UTC / Now / ±1 h / ±1 min | 관측 시각 입력·현재 시각·시각 이동 |
| Map | `model` 선택 |
| Coordinates | 적도/은하 좌표 표시 |
| Orbit / H / 북쪽 기준 위상 φ | 궤도·고도·위상 |
| Band / Model | 밴드와 Airglow 성분 `airglow_model` |
| NSIDE | 지도 해상도 |
| Min ELV | 최소 ELV |
| Evolve phase from epoch / Phase epoch UTC | 기준 시각부터 위상 자동 진행 |
| Sun/Moon avoidance / Min separation | 회피 적용과 각도 |
| Contour | contour 간격 |
| Calculate map | 현재 입력으로 지도 계산 |
| Save PDF / PNG / FITS.GZ | 현재 결과 저장 |

Map·Coordinates 선택과 시간 이동 버튼은 재계산을 요청합니다.
일반 입력칸을 수정한 뒤에는 **Calculate map**을 누릅니다.
Contour 간격 변경은 기존 결과를 다시 그립니다.
포인터를 지도 위로 옮기면 선택 좌표계의 좌표·밝기·천정각·ELV를 표시합니다.
계산은 작업 스레드에서 실행하며 첫 Numba 컴파일에는 시간이 더 걸릴 수 있습니다.
`--ra`, `--dec`는 CLI 좌표 질의 옵션이고 GUI에 표적 입력칸을 추가하지는 않습니다.

## 12. 빠른 계산과 반복 호출

`engine="auto"`는 Numba를 사용할 수 있으면 Airglow 병렬 적분을 사용하고,
없으면 NumPy 경로를 사용합니다. `engine="numpy"`로 명시할 수 있습니다.
`engine="numba"`인데 Numba가 없으면 오류를 냅니다.
Zodiacal/DGL만 선택하면 입력 지도를 내삽하며 Airglow LOS 적분은 수행하지 않습니다.

같은 `SkyModel`을 재사용하면 데이터 읽기와 해시 검증을 반복하지 않습니다.
처음에는 작은 NSIDE로 입력을 확인한 뒤 최종 결과에서 128로 높이는 방법이 편리합니다.
배열에 대한 반복 계산은 개별 `point()` 호출보다 `evaluate()`로 묶을 수 있습니다.

## 13. 좌표 보조 함수

| 함수 | 인터페이스 |
|---|---|
| `sdrift_sky.vector(ra, dec)` | degree 두 값을 적도 축 단위벡터 `(3,)`로 변환 |
| `sdrift_sky.coordinates(v)` | 단일 단위벡터를 `(경도, 위도)` degree로 변환 |
| `sky_coordinates.from_equatorial(v, coords)` | 적도 축 벡터 `(…,3)`를 선택 좌표계로 회전 |
| `sky_coordinates.to_equatorial(v, coords)` | 선택 좌표계 벡터를 적도 축으로 회전 |
| `sdrift_sky.utc_time(date)` | `(Unix 초, UTC 문자열)` 반환 |
| `sdrift_sky.ab_magnitude(flux)` | MJy/sr를 AB mag/arcsec²로 변환 |

`vector()`는 구면 경도·위도 벡터를 만드는 수식이므로 5.3절처럼 l, b를
변환할 때에도 사용할 수 있습니다. `coordinates()`는 벡터 자체의 축을 해석하지
않으므로 은하 축 벡터를 넣으면 l, b를 반환합니다.

## 14. 이전 코드에서 전환할 때

| 이전 방식 | 공통 API v2의 권장 방식 |
|---|---|
| `SkyModel`과 `AirglowOrbitModel` 사용법이 다름 | `from sdrift_sky import SkyModel`로 통일 |
| `kind="total"` | `model="total"` |
| `kind="airglow", model="glow"` | `model="airglow", airglow_model="glow"` |
| `--kind total --model combined` | `--model total --airglow-model combined` |
| Airglow의 `flux`만 광자 단위 | 모든 `flux`는 MJy/sr; 광자는 별도 필드 |
| 기존 SkyModel은 제외 방향에도 원래 밝기 반환 | 모든 모드에서 제외 방향은 NaN |
| 기존 SkyModel에 궤도 조건 없음 | 모든 모드에 기본 궤도·ELV 조건 적용 |
| 모드마다 FITS 열 구성 다름 | 공통 Python API v2는 같은 14개 열 |

`from airglow_orbit import AirglowOrbitModel`은 `SkyModel`의 호환 이름이며
**같은 새 API와 단위**를 사용합니다. 두 스크립트 CLI도 같습니다.
`kind=`, `phase=`, `--kind`, `--phase`와 구형 GLOW 모델 표기를 해석하는 호환
경로는 남겨 두었지만 새 코드에서는 권장 이름을 사용하십시오.
충돌하는 `model`과 `kind`를 동시에 주면 오류입니다.

이전 예제 `sky.point(0, 60, 날짜, band="VIS", model="zodiacal_dgl")`는 계속
실행됩니다. 다만 기본 궤도·ELV 조건과 NaN 정책은 새 규칙을 따릅니다.
과거 숫자와 비교할 때는 기본값·회피 조건·위상 정의·단위가 같은지 확인합니다.
이전 noon-midnight 위상 u가 정오 0° 기준이었다면 새 북쪽 기준 φ는 `(u−90°) mod 360°`입니다.

## 15. 자주 발생하는 문제

| 현상 | 확인할 내용 |
|---|---|
| `No module named ...` | 필요한 Python 파일과 requirements 설치 여부 |
| `manifest.json` 없음 | `docs/data/` 전체 또는 `data_dir` 경로 |
| SHA-256 불일치 | 코드와 데이터 버전 및 파일 손상; 데이터 파일을 함께 복원 |
| `Unknown sky model` 등 | `model`과 `airglow_model`의 역할, 대소문자 |
| 좌표 오류 | RA/DEC degree, DEC −90~90, 두 값을 함께 입력했는지 |
| 밝기가 NaN | `visible`, ELV, 태양·달 separation, 지구 차폐 확인 |
| 밝기가 무한대 | 계산 intensity가 0인지 확인 |
| 지도 전체가 비어 있음 | 해당 위상의 관측 영역과 ELV·회피각의 교집합 확인 |
| 빗금 표시 | Airglow 모델 지원 범위 또는 광전자 cutoff; 관측 금지의 뜻은 아님 |
| GUI가 열리지 않음 | 그래픽 데스크톱, Tkinter, Matplotlib backend |
| Numba 관련 오류 | Python/NumPy/Numba 호환 확인, 임시로 `engine="numpy"` 또는 `--engine numpy` |
| 예전 저장 파일과 열 번호가 다름 | `APIVER`, 열 이름과 단위를 확인 |

## 16. 물리 모델의 적용 범위

날짜는 천구상 위치, 태양·달 방향과 궤도 위상에 적용합니다.
Airglow 발광률은 평균 태양 활동의 사전 계산 표를 사용하며, 실시간 F10.7/Ap나
관측 시공간의 모든 대기 변화를 계산하는 예측 모델은 아닙니다.

DGL은 기존 C2022_all 상관 스펙트럼과 자기흡수 보정 가정을 사용하며
짧은 UV 파장으로의 연장은 관측값이 아닌 가정입니다.
Airglow에는 누락 성분, UV 흡수 미포함, dawn/dusk 화학 이력 미포함,
PALACE 야간층의 전 지구 적용 가정 등이 남아 있습니다.
달은 지구 중심 위치이며 위성 시차와 산란 달빛은 포함하지 않습니다.

추가 자료:

- [원본부터 데이터셋 생성](README_first.md)

- [궤도·ELV 기하와 배경 설명](README_airglow_orbit.md)
- [웹 실행 안내](https://github.com/seoncafe/SDRIFT_sky/blob/main/README.md)
- [Dawn–dusk 위상 그림](docs/assets/geometry_phase_dd.svg)
- [Noon–midnight 위상 그림](docs/assets/geometry_phase_nm.svg)
- [위성·지구 limb·ELV 기하 그림](docs/assets/geometry_limb_elv.svg)
- [GLOW 사용 조건](docs/data/glowlicense.txt)

# S-DRIFT Sky Background — web and Python

`docs/index.html`에서 zodiacal light, DGL, airglow를 계산합니다. 사용자는 Python,
healpy, NumPy를 설치할 필요가 없습니다. 최신 Chrome / Edge / Firefox / Safari에서
Web Worker, CompressionStream, DecompressionStream과 Web Crypto를 사용합니다.
외부 CDN이나 계산 서버 없이 작동하며, 데이터는 약 15 MB입니다.

## GitHub Pages

이 저장소의 **Settings → Pages → Build and deployment**에서
**Deploy from a branch → main (또는 실제 기본 브랜치) → /docs → Save**를 선택합니다.
배포 후 주소: https://seoncafe.github.io/SDRIFT_sky/

웹 실행에는 `docs/` 전체와 이 문서만 필요합니다. Python 사용자를 위한 `sdrift_sky.py`와 `requirements-python.txt`도 함께 제공합니다. 원본 FITS/NPZ,
Python GUI, GLOW/PALACE 소스, PDF 결과, 테스트와 생성 스크립트는 웹 실행에 필요하지
않습니다. `docs/.nojekyll`도 포함하십시오. 모든 자산 경로는 상대 경로입니다.

로컬 미리보기:

```sh
python3 -m http.server 8765 --directory docs
```

http://localhost:8765/ 에 접속하십시오. `file://`로 HTML을 직접 열면 데이터와
모듈 읽기가 차단될 수 있습니다. HTTPS 또는 localhost에서 사용하십시오.


## Python: RA, DEC, UTC 계산과 지도 표시

Python 3.10 이상에서 저장소를 내려받고 실행합니다. 원본 FITS나 GLOW/PALACE
설치 없이 저장소의 `docs/data/`를 사용합니다. 실행 중 인터넷 연결은 필요하지 않습니다.

```sh
python3 -m pip install -r requirements-python.txt
python3 sdrift_sky.py --ra 0 --dec 60 --date 2026-10-13T12:00:00Z --band VIS --model zodiacal_dgl --show
```

RA, DEC는 **degree**, 날짜는 **UTC**입니다. `--model`은 `zodiacal`, `dgl`,
`zodiacal_dgl`이며 기본값은 합산 모드입니다. `--band`는 UV1/UV2/UV3/VIS,
기본값은 VIS입니다. 날짜 범위는 1900–2100년입니다. 단일 좌표만 계산하려면
`--show`를 생략합니다. 표적은 지도에 파란 X로 표시됩니다.

화면 없이 파일 저장:

```sh
python3 sdrift_sky.py --ra 0 --dec 60 --date 2026-10-13T12:00:00Z \
  --band VIS --nside 128 --contour-step 0.2 \
  --pdf sky.pdf --png sky.png --fits sky.fits.gz --json point.json
```

세 모드 모두 Sun/Moon 위치, 빨간/주황색 점선 separation 경계, 회피 마스크를
공유합니다. 기본 Sun separation은 90°, Moon separation은 60°입니다.
`--sun-separation 100 --moon-separation 70`처럼 변경하거나
`--no-sun-avoidance`, `--no-moon-avoidance`로 각각 끕니다.
Contour는 기본 0.5 mag이고 `--contour-step`으로 0.1–1.0 mag를 선택합니다.
FITS.gz는 회피 영역을 NaN으로 저장합니다. 합산 FITS에는 두 성분도 저장됩니다.

좌표 질의는 주어진 방향에 대해 intensity를 내삽한 **회피 마스크 적용 전** 밝기와
실제 태양·달 separation, `visible` 여부를 반환합니다. 제외된 방향도 원래 밝기를
계산하며 지도에는 표시하지 않습니다. `display_ab_mag_arcsec2`는 이 경우 null입니다.
지도는 HEALPix 픽셀 중심에서 계산하므로 웹의 포인터 픽셀값과 좌표 직접 질의값은
픽셀 크기만큼 다른 방향을 나타낼 수 있습니다. 표시는 AB mag / arcsec²입니다.

Python에서 함수로 호출:

```python
from sdrift_sky import SkyModel
sky = SkyModel()
point = sky.point(0, 60, "2026-10-13T12:00:00Z", band="VIS", model="zodiacal_dgl")
print(point)
result = sky.sky_map("2026-10-13T12:00:00Z", nside=128, band="VIS", model="zodiacal_dgl")
figure = sky.plot(result, target=(0, 60), contour_step=0.5)
figure.savefig("sky.pdf", bbox_inches="tight")
```

`--show`에는 GUI backend와 그래픽 데스크톱이 필요합니다. 서버에서는
`--pdf`/`--png`를 사용합니다. Sun/Moon 방향은 웹과 동일한 저장된 지구 중심
GCRS 데이터이며 관측지·위성 시차나 달빛 산란은 포함하지 않습니다. DGL은
아래에 설명된 α ≈ α′ 보정 가정을 따르며 UV correlation 연장은 관측값이 아닙니다.

## 웹 사용

- 모델은 Zodiacal light / DGL / Zodiacal light + DGL / Airglow를 선택합니다.
- 모델과 UV1 (200–400 nm), UV2 (250–400 nm), UV3 (300–400 nm), VIS (400–700 nm)를 선택합니다.
- Sun avoidance를 켜면 Sun separation 경계를 빨간 점선으로 표시합니다. 태양 위치는 빨간 별과 Sun 라벨로 표시하며 Sun avoidance 설정을 꺼도 유지합니다.
- Zodiacal / DGL / Zodiacal + DGL 모두에 달의 위치와 회피각 경계를 표시합니다. Moon avoidance는 기본으로 활성화하며 최소 이각 60°, 5° 단위로 변경합니다. 태양과 달 회피를 동시에 켜면 두 조건을 모두 만족하는 방향만 표시·저장합니다. Airglow에는 날짜 기반 태양·달 마스크를 적용하지 않습니다.
- 달은 Astropy builtin으로 계산한 지구 중심 GCRS 위치를 3시간 간격으로 저장하고 내삽합니다. 1900–2100년을 지원하며 관측지/위성 시차 및 달빛 산란은 포함하지 않습니다. FITS에는 MOONMASK, MOONMIN, MOONRA, MOONDEC, MOONOBS를 기록합니다.
- 기본 모델은 Zodiacal light + DGL이며 시작 시 현재 UTC로 계산합니다. 날짜·시간 옆 ▲/▼와 키보드 ↑/↓는
  월·연도, 일·월, 초·분 등의 경계를 넘겨 변경합니다. Sun separation은 기본 90°, 5° 단위입니다.
- DGL과 Zodiacal light + DGL은 zodiacal과 같은 날짜·시간 및 Sun avoidance 설정을
  사용합니다. DGL 분포는 RA, DEC에 고정되며 두 성분의 intensity를 더해 AB magnitude로 변환합니다.
  합산 FITS에는 각 성분 intensity도 포함됩니다.
- Airglow는 H와 위성 위치 SZA를 지정합니다. 기본 H=550 km, SZA=90°,
  GLOW + PALACE night layers입니다. 평균 태양 활동 모델로 날짜는 사용하지 않습니다.
- NSIDE 기본값은 128입니다. 저성능 기기에서는 32 또는 64로 선택하십시오.
  변경 시 자동 계산하며, `중지` 버튼으로 취소할 수 있습니다.
- RA는 왼쪽으로 증가하며 모든 지도에서 contour 간격은 기본 0.5 mag이며, 0.1–1.0 mag를 0.1 mag 단위로 선택할 수 있습니다. 선택을 바꾸면 기존 지도에 바로 적용되고 PDF/PNG도 같은 간격으로 저장됩니다.
  큰 AB magnitude의 색은 어둡습니다. 지도 위 포인터로 RA, DEC와 값을 읽습니다.
- Contour level마다 기존 방식으로 라벨을 표시합니다. 낮은 level에서 높은 level 순으로 선굵기를 0.6–2.2 px로 증가시키며, 화면과 PDF/PNG 저장에 동일하게 적용됩니다.
- FITS.gz: RING / Equatorial HEALPix BINTABLE. Airglow에는 EARTH_OCCULTED,
  MODEL_VALID, MIN_ALTITUDE, PHOTO_CUTOFF와 sunlit radiance도 저장합니다.
- PDF 저장 / 인쇄: 인쇄 창에서 대상을 PDF로 선택하십시오. 지도 제목·범례와
  모델 라이선스 표기가 포함됩니다. PNG도 저장할 수 있습니다.
- Zodiacal 파일명은 `zodiacal_VIS_nside128_20261008_123000_anti-sun.fits.gz`
  형식입니다. 시간은 UTC이며 anti-sun은 마지막에 붙습니다.

## 계산의 의미

Constant **photon-counting** throughput에 대한 AB mag/arcsec²입니다.
Zodiacal intensity를 회전·내삽한 뒤 AB magnitude로 변환하며,
NSIDE 128 춘분 기준 template을 사용합니다. 날짜별 태양 벡터는 Python과 같은
Astropy get_sun의 J2000 축 벡터를 12시간 간격으로 저장해 내삽합니다.
지원 날짜는 UTC 1900–2100년입니다. Sun avoidance는 Sun separation ≥ 지정값인 방향만 남기며 기본값은 90°입니다.

DGL은 IRIS 100 μm와 SFD E(B−V), Fitzpatrick 1999 R_V=3.1을 사용합니다.
C2022_all 표의 Å를 μm로 변환하여 보간하고 원본 보정계수 2.1을 유지합니다.
372.5 nm 미만의 correlation은 0.127596으로 일정하게 연장합니다.
논문 식 (3)–(6)의 자기흡수식 `I_nu = alpha * beta * I100 * lambda_um/100`을 적용합니다.
단, C2022_all은 Figure 3의 alpha_prime 표이며 원래 alpha 및 표본 유효 beta를
확보하지 못했으므로, **beta_ref=1 (광학적으로 얇은 표본), alpha ≈ alpha_prime**을
가정합니다. 이는 자기흡수 법칙은 논문에 따른 계산이지만 정규화는 근사입니다.
보정계수 2.1은 유지하며, 기존 inverse-beta 지도는 사용하지 않습니다.
UV 연장 및 은하면/높은 광학 깊이에 대한 외삽은 가정입니다.
파장 적분은 constant photon-counting, 최대 0.5 nm 간격이며,
DGL은 날짜에 따라 회전하지 않습니다. 색상 범위는 마스크 적용 후 유효한 mag의 최솟값·최댓값으로 자동 설정하고,
0.5 mag 경계로 바깥쪽 반올림합니다.
FITS에 BETAOP=MULTIPLY, ALPHAPP=T, BETAREF=1을 기록합니다.
출처: https://arxiv.org/abs/2201.01378 .

Airglow는 기존 q(h, SZA, component)를 band별 day/all 및 night-only 두 채널로
합산한 데이터를 사용합니다. 흡수가 없는 constant throughput에서는 이 합산과
선별 합산 적분이 동등합니다. 브라우저에서 고도·국소 SZA를 쌍선형 내삽하고,
각 HEALPix 시선에 2 km midpoint 적분을 합니다. Python의 H/SZA 적분 LUT를
재내삽하지 않으므로 그 LUT의 추가 오차는 없습니다. Float32 저장 반올림은 있습니다.

태양은 RA=0°, DEC=0°, 위성 위치는
`(R+H) × (cos(SZA), 0, sin(SZA))`로 고정합니다. 실제 날짜의 위성 궤도나
자세에 연결된 관측 RA/DEC가 아니라 이 기준 geometry의 equatorial map입니다.
R=6371 km, 데이터 상한은 1500 km입니다. 태양 그림자에서만 PALACE 밤 성분을
사용합니다. 시간에 따른 화학 반응과 twilight history는 적용하지 않습니다.

Airglow band 성분은 **불완전**합니다. 특히 UV1/UV2의 200–300 nm 분자·연속
스펙트럼과 UV LOS 흡수, H geocorona, 일부 고진동 준위/회전 스펙트럼이
없습니다. PALACE는 Paranal 밤 기둥값을 전 지구 그림자 안에 적용하고 고도
분포를 가정합니다. F10.7=100, Ap=4, 적도 4계절×4 LST 평균이며 전 지구
관측 climatology가 아닙니다. 회색은 지구 차폐, 빗금은 최소 지원 고도 또는
GLOW photoelectron cutoff(SZA ≥ 1.85 rad)로 인해 모델 적용이 불완전한 방향입니다.
0 radiance는 +∞ mag이고 유한 색 범위/contour에서는 제외합니다.

출처·성분·입력값·SHA-256은 `docs/data/manifest.json`에 포함되어 있습니다.
PALACE: https://gmd.copernicus.org/articles/18/4353/2025/ (CC BY 4.0).
GLOW와 zodiacal 라이선스는 `docs/data/`에 포함되어 있습니다.

This software is based on calculations performed using the GLOW model.
Use of the GLOW model is governed by the Open Source Academic Research License
Agreement contained in [glowlicense.txt](docs/data/glowlicense.txt).
Browser routines modified 2026-10-08.


## 두 이상적 궤도의 airglow 지도

Airglow 모드에서 dawn-dusk / noon-midnight, UTC, 고도, 위상, 최소 ELV를
선택합니다. 두 궤도 모두 태양·달 위치와 경계를 표시하며 기본 Sun/Moon separation은
90°/60°입니다. Python GUI는 `python3 airglow_orbit.py --gui`로 실행합니다.
FITS.gz, PDF, PNG 저장과 CLI 계산도 지원합니다.
자세한 기하, 위상 자동 진행, 실행 예제와 제한은
[README_airglow_orbit.md](README_airglow_orbit.md)를 참조하십시오.

Zodiacal light + DGL + Airglow 합산 모드는 세 성분의 intensity를 더하며
Airglow와 동일한 궤도·각도·ELV 및 태양·달 회피 조건을 사용합니다.
웹 모델 목록에서 선택하거나 `python3 airglow_orbit.py --gui --kind total`로
실행하십시오. [사용법](README_airglow_orbit.md#zodiacal-light--dgl--airglow-합산).

웹과 Python 궤도 뷰어의 기본 모드는 **Zodiacal light + DGL + Airglow**입니다.

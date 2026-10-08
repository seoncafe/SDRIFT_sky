# S-DRIFT Sky Background — browser app

`docs/index.html`에서 zodiacal light, DGL, airglow를 계산합니다. 사용자는 Python,
healpy, NumPy를 설치할 필요가 없습니다. 최신 Chrome / Edge / Firefox / Safari에서
Web Worker, CompressionStream, DecompressionStream과 Web Crypto를 사용합니다.
외부 CDN이나 계산 서버 없이 작동하며, 데이터는 약 8 MB입니다.

## GitHub Pages

이 저장소의 **Settings → Pages → Build and deployment**에서
**Deploy from a branch → main (또는 실제 기본 브랜치) → /docs → Save**를 선택합니다.
배포 후 주소: https://seoncafe.github.io/SDRIFT_sky/

GitHub Pages의 공식 설정 안내:
https://docs.github.com/en/pages/getting-started-with-github-pages/creating-a-github-pages-site

사이트에 필요한 파일은 `docs/` 전체와 이 문서뿐입니다. 원본 FITS/NPZ,
Python GUI, GLOW/PALACE 소스, PDF 결과, 테스트와 생성 스크립트는 웹 실행에 필요하지
않습니다. `docs/.nojekyll`도 포함하십시오. 모든 자산 경로는 상대 경로입니다.

로컬 미리보기:

```sh
python3 -m http.server 8765 --directory docs
```

http://localhost:8765/ 에 접속하십시오. `file://`로 HTML을 직접 열면 데이터와
모듈 읽기가 차단될 수 있습니다. HTTPS 또는 localhost에서 사용하십시오.

## 사용

- 모델은 Zodiacal light / DGL / Zodiacal light + DGL / Airglow를 선택합니다.
- 모델과 UV1 (200–400 nm), UV2 (250–400 nm), UV3 (300–400 nm), VIS (400–700 nm)를 선택합니다.
- Zodiacal light는 시작 시 현재 UTC로 계산합니다. 날짜·시간 옆 ▲/▼와 키보드 ↑/↓는
  월·연도, 일·월, 초·분 등의 경계를 넘겨 변경합니다. 태양 이각은 기본 60°, 5° 단위입니다.
- DGL과 Zodiacal light + DGL은 zodiacal과 같은 날짜·시간 및 anti-sun 설정을
  사용합니다. DGL 분포는 RA, DEC에 고정되며 두 성분의 intensity를 더해 AB magnitude로 변환합니다.
  합산 FITS에는 각 성분 intensity도 포함됩니다.
- Airglow는 H와 위성 위치 SZA를 지정합니다. 기본 H=550 km, SZA=90°,
  GLOW + PALACE night layers입니다. 평균 태양 활동 모델로 날짜는 사용하지 않습니다.
- NSIDE 기본값은 128입니다. 저성능 기기에서는 32 또는 64로 선택하십시오.
  변경 시 자동 계산하며, `중지` 버튼으로 취소할 수 있습니다.
- RA는 왼쪽으로 증가하며 모든 지도에서 contour 간격은 0.5 mag입니다.
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
지원 날짜는 UTC 1900–2100년입니다. Anti-sun은 태양과의 시선각 ≥ 지정값입니다.

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

This software is part of the GLOW model. Use is governed by the Open Source
Academic Research License Agreement contained in the file
[glowlicense.txt](docs/data/glowlicense.txt).
Browser routines modified 2026-10-08.


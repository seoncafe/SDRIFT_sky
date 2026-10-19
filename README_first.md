# S-DRIFT 데이터셋을 처음부터 만드는 방법

작성 기준: 2026-10-19. 아래 명령은 **생성 스크립트가 있는 개발 폴더의 루트**에서 실행합니다.
UV1, UV2, UV3, VIS의 기준 데이터를 만들고, 현재 웹 및 공통 Python API가 읽는
`docs/data/`로 변환하는 순서입니다. 기존 데이터를 지우거나 계산을 시작하는 문서가
아니며, 필요한 단계의 명령을 직접 실행합니다.

**GitHub 사용자 배포본만으로는 원본부터 재생성할 수 없습니다.**
[SDRIFT_sky 저장소](https://github.com/seoncafe/SDRIFT_sky)는 계산 완료 데이터와
조회 프로그램을 배포하며, 아래 생성 스크립트·GLOW/PALACE 원본·대용량 IRIS/SFD
자료는 포함하지 않습니다. 처음부터 만들려면 전체 개발 폴더의 백업과 원본 자료를
준비해야 합니다. 이미 배포된 데이터로 지도만 그리려면
[Python 인터페이스 안내서](README_python_interface.md)를 사용하십시오.

## 1. 전체 생성 순서와 밴드

```text
skysurf 원본 ── zodiacal_healpix.py ── zodi_data/*_template.fits.gz ──┐
IRIS + SFD + C2022_all ── dgl_healpix.py ── dgl_data/*_template.fits.gz ┤
GLOW ── prepare_airglow_glow.py ── glow_selected_lines_mean.npz ──────┤
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
                                             공통 Python API / GUI / Web
```

Airglow는 날짜 하나의 기준 전천 지도를 회전하는 방식이 아닙니다.
고도 `h`, 방출 지점의 태양천정각 `chi`, 발광 성분에 따른 방출률
`q(chi, h, component)`를 먼저 저장합니다. 관측 날짜·궤도·위성 고도·위상을
지정하면 각 시선에서 이 표를 내삽하고 적분하여 지도를 만듭니다.

| 밴드 | 파장 범위 |
|---|---|
| UV1 | 200–400 nm |
| UV2 | 250–400 nm |
| UV3 | 300–400 nm |
| VIS | 400–700 nm |

모두 constant throughput 및 photon-counting AB 가정입니다.
Zodiacal/DGL 기준 지도는 RING HEALPix, NSIDE=128, 적도 좌표이며 픽셀은
196,608개입니다. 기준 지도에는 태양·달 회피 마스크를 넣지 않습니다.
관측 시점의 회피 조건과 최소 ELV는 최종 지도 계산 때 적용합니다.

## 2. 개발 코드와 실행 환경 준비

### 2.1 필요한 개발 파일

다음 파일을 전체 개발 폴더 백업에서 복원합니다. 원본 모델만 내려받아서는
S-DRIFT의 생성·변환 작업을 수행할 수 없습니다.

```text
SDRIFT_sky/
├── zodiacal_healpix.py
├── dgl_healpix.py
├── dgl_model/                    # 수정된 DGL_lib.py, fitzpatrick.py, __init__.py
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
├── docs/assets/                  # 웹 코드
├── docs/index.html
├── skysurf/                      # 3절에서 준비
├── GLOW/                         # 5절에서 준비
├── PALACE/                       # 6절에서 준비
├── dgl_data/                     # 4절의 원본 자료
├── zodi_data/                    # 생성 결과
└── airglow_data/                 # 생성 결과 및 빌드 로그
```

선택적인 기존 H/SZA 지도와 LOS LUT 생성에는 `airglow_healpix.py`,
`prepare_airglow_lut.py`, `airglow_lut.py`도 필요합니다.
검증 스크립트는 10절에 나열했습니다.

### 2.2 Python 패키지와 컴파일러

별도 환경을 만들고 같은 Python으로 설치·생성을 수행합니다.
`python3`는 그 환경의 Python을 가리켜야 합니다.

```bash
python3 -m venv .venv-build
source .venv-build/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install "numpy>=2.0" scipy astropy healpy matplotlib Cython setuptools
mkdir -p zodi_data dgl_data airglow_data
```

Python 3.10 이상을 사용합니다. 생성 코드에서 `np.trapezoid()`를 사용하므로
여기서는 **NumPy 2.0 이상**을 명시했습니다. 오래된 requirements의 하한만
만족하는 환경과 구분하십시오. Cython은 PALACE의 순수 Python 실행에도 필요합니다.
`zodiacal_healpix.py`는 모델 폴더를 직접 import하므로 별도 zodisurf 설치는 필요하지 않습니다.

GLOW 테이블 생성에는 PATH에서 실행 가능한 `gfortran`이 필요합니다.
macOS Homebrew를 사용하는 경우의 설치 예입니다.

```bash
brew install gcc
gfortran --version
```

Numba와 MPI는 필수가 아닙니다. 우선 NumPy 및 로컬 프로세스 방식으로 생성할 수 있습니다.
Numba를 쓰려면 현재 Python·NumPy와 호환되는 버전을 설치하고 import를 확인합니다.
특정 과거 환경의 버전 충돌 여부를 현재 환경에도 그대로 적용하지 마십시오.

```bash
python3 -m pip install numba
python3 -c "import numpy, numba; print(numpy.__version__, numba.__version__)"
```

GUI에는 Python 버전에 맞는 Tkinter/Tcl/Tk와 그래픽 데스크톱이 추가로 필요합니다.
파일 생성만 하는 단계에는 Tkinter가 필요하지 않습니다.

## 3. Zodiacal light: 춘분 기준 FITS 생성

### 3.1 원본 모델 준비

원본은 [Zodiacal light 모델 저장소](https://github.com/rosaliaobrien/skysurf)의
`Zodiacal_Light_Model_Python/`입니다. 이 프로젝트의 로컬 디렉터리 이름은
호환성을 위해 `skysurf/`를 유지합니다.

새 폴더에서만 아래 명령으로 받습니다. 이미 원본이 있으면 다시 clone할 필요가 없습니다.

```bash
git clone https://github.com/rosaliaobrien/skysurf.git skysurf
```

현재 개발 폴더에서 확인한 원본 commit은
`92dad09e18b9b69e51035a0fdf51364cb67e6ea6`입니다.
동일한 원본 버전이 필요하면 새 clone에서 다음과 같이 고정합니다.

```bash
git -C skysurf checkout 92dad09e18b9b69e51035a0fdf51364cb67e6ea6
```

`skysurf/Zodiacal_Light_Model_Python/zodisurf/get_zmod.py`와 `LICENSE`가 있어야 합니다.

### 3.2 네 밴드 생성

```bash
python3 zodiacal_healpix.py --make-template --nside 128 --band all \
  --reference-date 2000-03-20T07:35:00 --workers 4 --no-plot
```

`--workers`는 로컬 프로세스 수입니다. 자원에 맞게 줄이거나 늘립니다.
이 단계는 태양·달 마스크가 없는 전천 FITS를 생성합니다.
`--no-plot`은 불필요한 PDF 생성을 생략합니다.

```text
zodi_data/zodiacal_UV1_nside128_template.fits.gz
zodi_data/zodiacal_UV2_nside128_template.fits.gz
zodi_data/zodiacal_UV3_nside128_template.fits.gz
zodi_data/zodiacal_VIS_nside128_template.fits.gz
```

각 파일의 두 열은 `AB_MAG_ARCSEC2`, `INTENSITY_MJY_SR`입니다.
기준 UTC는 2000년 춘분인 `2000-03-20 07:35:00 UTC`이며,
밴드 내 파장 간격은 0.01 μm입니다. 날짜별 회전은 태양의 황경 변화로 계산합니다.
UV1의 200–250 nm 구간에서는 원본 위상함수의 250 nm 매개변수를 사용하는
제한이 있습니다.

한 밴드만 만들려면 `--band VIS`처럼 지정합니다.
이 스크립트는 기존 기준 파일을 덮어쓰므로 보존할 결과는 실행 전에 별도 보관합니다.
웹 변환기는 NSIDE=128의 위 파일명 네 개를 읽습니다.

## 4. DGL: C2022_all 기준 FITS 생성

### 4.1 입력 파일 준비

다음 자료를 `dgl_data/`에 둡니다. 현재 개발 폴더에는 모두 복사되어 있습니다.
`../DGL` 등을 가리키는 symlink는 필요하지 않습니다.

| 파일 | 내용 |
|---|---|
| `DGL_C2022_all.txt` | C2022 Fig. 3 오른쪽 correlation 표; 파장 단위 Å |
| `IRIS_nohole_4_2048_v2.fits.gz` | 100 μm IRIS, Galactic RING, NSIDE=2048 |
| `lambda_sfd_ebv.fits.gz` | SFD E(B−V), Galactic NESTED, NSIDE=512 |
| `SOURCE_README.txt` | correlation 표의 원본 설명 |
| `source_files.json` | 입력 자료의 출처·복사 경로·SHA-256 기록 |

현재 자료의 원본 위치는 `dgl_data/source_files.json`에 기록되어 있습니다.
기존 작업 환경에서는 표와 FITS를 `../K-DRIFT_sim/DGL_data/`에서 복사했습니다.
새 환경에서는 **백업의 동일 파일**을 복원하는 것이 현재 결과 재현 경로입니다.
IRIS라는 이름의 다른 지도나 논문에서 새로 digitize한 표가 자동으로 동일 입력이
되는 것은 아닙니다. 특히 이 IRIS 가공 파일과 C2022_all 표에 대한 직접 다운로드
절차는 현재 생성 코드에 구현되어 있지 않습니다.

`dgl_model/`은 프로젝트에서 수정한 복사본을 유지합니다. 원래 `../DGL/DGL_lib.py`로
덮어쓰면 파장 단위 보정과 논문 방식의 자기흡수 처방이 바뀔 수 있습니다.
`source_files.json`의 코드 hash는 복사 전 원본 기록이므로 현재 수정본 hash와
같아야 하는 것은 아닙니다.

자료를 새로 교체했다면 `source_files.json`도 실제 입력에 맞게 갱신합니다.
원본 필드는 기록이며 계산 중 그 절대 경로에 접근하지는 않습니다.

### 4.2 생성 명령과 물리 가정

```bash
python3 dgl_healpix.py --nside 128 --step-nm 0.5
```

기존 파일을 다시 만들 때만 `--overwrite`를 추가합니다.

```text
dgl_data/dgl_UV1_nside128_template.fits.gz
dgl_data/dgl_UV2_nside128_template.fits.gz
dgl_data/dgl_UV3_nside128_template.fits.gz
dgl_data/dgl_VIS_nside128_template.fits.gz
dgl_data/dgl_nside128_metadata.json
```

출력은 적도 좌표의 RING HEALPix이며 두 열의 단위는 Zodiacal과 같습니다.
원본 FIR/E(B−V)를 먼저 NSIDE=128로 낮추고 Galactic에서 계산한 intensity를
적도 좌표 픽셀 중심으로 내삽합니다. 날짜에 따라 DGL 자체는 회전하지 않습니다.

현재 처방은 다음과 같습니다.

```text
alpha_approx(lambda) = 2.1 × C2022_all(lambda)
beta = (1 − exp(−tau)) / tau
I_nu(lambda) = alpha_approx × beta × I100 × lambda_um / 100
<I_nu> = integral[I_nu(lambda) d(lambda)/lambda] / ln(lambda_max/lambda_min)
```

[Chellew et al. (2022)](https://arxiv.org/abs/2201.01378)의 자기흡수식처럼
**beta를 곱합니다.** 원본의 inverse-beta 처방으로 돌아가지 않습니다.
다만 표는 alpha_prime이므로 `alpha ≈ alpha_prime`, `beta_ref=1`이라는
정규화 근사를 사용합니다. 원 논문 표본의 intrinsic alpha를 정확히 복원한 것은 아닙니다.
372.5 nm 미만은 가장 짧은 파장의 correlation 값을 일정하게 연장하되,
Fitzpatrick 소광과 파장에 따른 intensity 변환은 유지합니다.

## 5. Airglow 1단계: GLOW 방출률 표 생성

### 5.1 GLOW 원본 준비

원본 안내는 [NCAR/HAO GLOW](https://www2.hao.ucar.edu/modeling/glow)에 있습니다.
현재 생성기는 로컬 **GLOW v0.982**의 Fortran 파일 구조를 기준으로 작성했습니다.
동일 결과 재현에는 개발 폴더에 보관한 이 원본과 `data/`를 복원합니다.
새 GLOW 버전을 받는 경우에는 빌드 패치가 호환되는지 확인해야 합니다.

```text
GLOW/
├── Makefile.glowbasic.gfort
├── Glowlicense.txt
├── Glow.txt
├── cglow.f90, glow.f90, mzgrid.f90, gchem.f90, ...
└── data/
```

공식 안내 페이지는 원본 입수 및 라이선스 확인 출처입니다. 특정 최신 archive
파일명을 가정한 자동 다운로드 명령은 제공하지 않습니다.
기존 `GLOW/glow.exe`를 직접 실행하는 단계도 필요하지 않습니다.

### 5.2 평균 입력 테이블 생성

```bash
python3 prepare_airglow_glow.py --f107 100 --ap 4 --chi-step 1 \
  --top-km 1500 --radius-km 6371
```

생성 결과:

```text
airglow_data/glow_selected_lines_mean.npz
airglow_data/glow_selected_lines_mean.json
airglow_data/build_glow/compile.log
airglow_data/build_glow/table.log
```

NPZ의 주요 키는 `altitude_km`, `chi_deg`, `wavelength_nm`, `q`, `metadata`입니다.
`q`의 배열 순서는 `(Nchi, Nz, Ncomponent)`, 단위는 photons/cm³/s입니다.
`q`는 4π로 나누기 전의 체적 광자 생성률이며, 시선 적분에서 `1/(4π)`를 적용합니다.
`component`에는 개별 원자선과 대표 파장으로 처리하는 분자 채널 등이 포함됩니다.

기본 조건은 F10.7=100 sfu, Ap=4, 오로라 입력 0, GLOW xuvfac=3이며,
SZA 0–180°를 1° 간격으로 계산합니다. GLOW 기본 고도 격자를 1500 km까지
확장합니다. 배경 대기는 적도·경도 0°의 네 계절 및 네 지방시 입력장을
평균한 기준 대기입니다. **평균 입력으로 계산한 결과이며 전 지구 기후 평균은 아닙니다.**

생성기는 `airglow_data/build_glow/`에서 파생 소스를 컴파일하고 원본을 보존합니다.
이 빌드 안의 `data`는 원본 `GLOW/data`를 가리키므로 원본 폴더를 옮기면 다시
빌드해야 합니다. 기본 NPZ/JSON은 재실행 시 덮어씁니다.

## 6. Airglow 2단계: PALACE 야간 성분 결합

### 6.1 원본 archive와 디렉터리

PALACE v1.0 원본 배포 식별자는
[Zenodo DOI 10.5281/zenodo.14064022](https://doi.org/10.5281/zenodo.14064022)입니다.
현재 개발 폴더에는 아래 세 archive가 보관되어 있습니다.

```text
PALACE/PALACE.zip
PALACE/test.zip
PALACE/PMD.zip
```

모델 계산에 쓰는 것은 `PALACE.zip`의 코드·데이터이며, `test.zip`은 배포 예제
비교용입니다. `PMD.zip`은 보조 자료이지만 현재 `verify_palace.py`가 세 archive
모두의 CRC를 검사하므로 **그 검증을 실행하려면 세 개가 모두 필요합니다.**

아직 풀지 않은 새 환경에서만 다음 예처럼 풉니다. 이미 `model/`,
`reference_test/`가 있으면 그대로 사용합니다.

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

최종적으로 다음 경로에 파일이 있어야 합니다.

```text
PALACE/model/src/palace/palace.py
PALACE/model/src/palace/config/palace_default.par
PALACE/model/src/palace/data/palace_lines.fits
PALACE/model/src/palace/data/palace_cont.fits
PALACE/model/src/palace/data/palace_var.fits
PALACE/reference_test/palace_test.par
PALACE/reference_test/palace_test.fits
```

`prepare_airglow_palace.py`는 위 로컬 모델을 직접 import합니다.
가상환경에 Cython이 있으면 `setup_cython.py`로 확장 모듈을 빌드할 필요는 없습니다.
이전 `PALACE/dependencies/`의 Cython 복사본을 새 환경에 가져올 경우에는
가상환경의 패키지보다 import 우선순위가 높다는 점을 확인합니다.

### 6.2 배포 예제 검증과 결합

```bash
python3 verify_palace.py
python3 prepare_airglow_palace.py
```

결합 테이블이 이미 있으면 마지막 명령에 `--overwrite`를 추가합니다.
검증은 `PALACE/verification.json`과
`airglow_data/palace_night_reference_300_700nm.fits.gz`를 만듭니다.
참조 스펙트럼은 지상 천정의 야간 평균이며 위성 전천 지도가 아닙니다.

결합 출력은 다음과 같습니다.

```text
airglow_data/glow_palace_night_mean.npz
airglow_data/glow_palace_night_mean.json
airglow_data/glow_palace_night_mean_atomic_manifest.json
airglow_data/glow_palace_night_mean_molecular_manifest.json
```

PALACE의 OH/O2/Na, 연속광, 중복 없는 O선, 중간권 557.7 nm 잔차와 297.2 nm
분기를 추가하며 N2/N2+의 같은 상부 진동 준위에서 나오는 분자 밴드도 결합합니다.
분자 전이계수는 `airglow_molecular.py`에 있으므로 NIST PDF를 실행 중 읽지 않습니다.
PALACE의 기둥 세기는 가정한 Gaussian 고도 분포로 바꾸고, 추가 야간 성분은
각 방출 지점이 기하학적 지구 그림자에 들어가는 경우에 적용합니다.

기본 FWHM은 OH/녹색선 8 km, 기타 중간권 10 km, 열권 원자선 100 km입니다.
완전한 관측 고도 분포의 복원이 아닙니다. PALACE의 300 nm 미만, H 지구코로나,
완전한 N2 LBH, UV 흡수, 황혼 화학 이력 등은 미구현입니다.
최신 성분 판정은 로컬 `AIRGLOW_COMPONENTS.md`를 참조하십시오.

## 7. 선택 단계: Airglow FITS와 LOS LUT

**현재 공통 Python API와 웹 실행을 위한 필수 단계는 아닙니다.**
8절의 웹 변환기는 두 방출률 NPZ를 읽으며, 아래 지도나 `*_los_lut.npz`는 읽지 않습니다.
현재 궤도 뷰어는 밴드별 방출률 표를 내삽하고 LOS 적분합니다.

### 7.1 H/SZA 고정 기준 지도

이전 기준 기하에서 결과를 점검하거나 FITS로 보관할 때 사용합니다.

```bash
python3 airglow_healpix.py --table airglow_data/glow_selected_lines_mean.npz \
  --H 550 --sza 90 --band UV1 UV2 UV3 VIS --nside 128 --workers 4 --no-pdf
python3 airglow_healpix.py --table airglow_data/glow_palace_night_mean.npz \
  --H 550 --sza 90 --band UV1 UV2 UV3 VIS --nside 128 --workers 4 --no-pdf
```

예를 들어 VIS 출력은 각각 다음과 같습니다.

```text
airglow_data/airglow_VIS_nside128_template_H550_sza90.fits.gz
airglow_data/airglow_VIS_nside128_template_H550_sza90_glow_palace.fits.gz
```

이 지도의 기하는 H와 SZA를 고정한 기준 기하이며 특정 날짜의 궤도 위치를 뜻하지
않습니다. 위 명령의 H=550 km는 점검 사례이고 현재 궤도 뷰어 기본 고도는 650 km입니다.
기존 파일 재계산에는 `--overwrite`를 추가합니다.

### 7.2 이전 H/SZA 뷰어용 적분 LUT

```bash
python3 prepare_airglow_lut.py --table airglow_data/glow_selected_lines_mean.npz --workers 4
python3 prepare_airglow_lut.py --table airglow_data/glow_palace_night_mean.npz --workers 4
```

각각 `glow_selected_lines_mean_los_lut.npz`,
`glow_palace_night_mean_los_lut.npz`가 같은 폴더에 생성됩니다.
많은 H/SZA/방위각/접선 고도 조합의 적분이므로 이 단계는 오래 걸릴 수 있습니다.
방출률 표를 바꾸면 해당 LUT도 다시 만들어야 합니다. 원본 SHA-256이 다르면
기존 LUT를 그대로 사용할 수 없습니다. 재생성에는 `--overwrite`를 추가합니다.

선택한 LUT의 오차는 다음 명령으로 직접 적분과 비교합니다.

```bash
python3 validate_airglow_lut.py --table airglow_data/glow_palace_night_mean.npz --nside 16
```

8절로 바로 넘어가도 현재 공통 Python API와 웹 지도를 사용할 수 있습니다.

## 8. Python·웹 공용 데이터로 변환

다음 자료가 모두 준비된 뒤 실행합니다.

- `zodi_data/zodiacal_{UV1,UV2,UV3,VIS}_nside128_template.fits.gz`
- `dgl_data/dgl_{UV1,UV2,UV3,VIS}_nside128_template.fits.gz`
- `dgl_data/dgl_nside128_metadata.json`, correlation 표와 원본 설명
- `airglow_data/glow_selected_lines_mean.npz`
- `airglow_data/glow_palace_night_mean.npz`
- `GLOW/Glowlicense.txt`, `skysurf/Zodiacal_Light_Model_Python/LICENSE`

```bash
python3 export_web_data.py
```

이 단계는 다음을 생성합니다.

| 출력 | 용도 |
|---|---|
| `docs/data/manifest.json` | 밴드·차원·격자·출처·SHA-256·모델 가정 |
| `zodiacal_{band}.f32.gz` | Zodiacal의 MJy/sr 기준 지도 |
| `dgl_{band}.f32.gz` | DGL의 MJy/sr 기준 지도 |
| `airglow_{glow,combined}_{band}.f32.gz` | `(Nchi, Nz, 2)` 밴드 방출률 표 |
| `sun_j2000.f32.gz` | 12시간 간격 태양 방향 |
| `moon_gcrs.f32.gz` | 3시간 간격 달 방향 |
| 라이선스·출처 파일 | GLOW, Zodiacal, DGL 자료 설명 |

표의 파일명은 모두 `docs/data/` 아래이며 `{band}`는 네 밴드 각각을 뜻합니다.
Float32 little-endian 배열을 gzip한 파일입니다. `.f32.gz`는 FITS가 아닙니다.
Airglow 마지막 축의 2개 값은 일반 성분과 그림자 안에서만 적용할 야간 성분의
밴드별 광자 생성률이며, 아직 시선 적분을 한 radiance가 아닙니다.
태양·달 표는 Astropy builtin으로 1900–2100년 사용 범위를 생성합니다.
월력 전체를 계산하므로 변환도 시간이 걸릴 수 있습니다.

**DGL 입력이 빠져도 exporter가 DGL을 생략하고 끝날 수 있습니다.**
전체 모델을 배포하려면 DGL을 먼저 생성하고, 10절에서 manifest에 `dgl`이 있는지
반드시 확인합니다. 변환기는 NSIDE=128의 고정 파일명을 사용하므로 다른 해상도
파일만 만들어 놓으면 현재 명령으로는 변환할 수 없습니다.

DGL만 바꾼 경우 기존 manifest를 유지하며 갱신할 수도 있습니다.

```bash
python3 export_web_data.py --dgl-only
```

이 옵션과 `--moon-only`는 기존 manifest가 있어야 합니다. 원본부터 만드는 첫 실행은
옵션 없는 전체 변환을 사용합니다. Zodiacal 또는 Airglow를 바꾼 경우에도 전체 변환을
사용합니다. exporter는 출력 파일과 manifest를 덮어쓰므로 서비스 중인 폴더보다는
개발 폴더에서 생성·검증한 뒤 배포합니다.

## 9. 최종 지도와 합산 모드 확인

현재 공통 인터페이스를 사용합니다. 기준 자료 생성 스크립트의 예전 CLI 옵션과
혼동하지 않도록 최종 관측 지도의 모드는 `--model`, Airglow 내부 성분은
`--airglow-model`로 선택합니다.

```bash
python3 sdrift_sky.py --model total --airglow-model combined --band VIS \
  --date 2026-10-19T00:00:00Z --orbit dawn-dusk --height 650 \
  --orbital-phase-deg 0 --min-elv 20 --sun-separation 90 --moon-separation 60 \
  --nside 128 --contour-step 0.2 \
  --fits output/total_VIS_nside128_20261019.fits.gz \
  --pdf output/total_VIS_nside128_20261019.pdf
```

`total`은 Zodiacal + DGL + Airglow입니다. 세 성분을 같은 밴드의 **MJy/sr로 더한 뒤**
AB mag/arcsec²로 변환하며 등급끼리 더하지 않습니다.
최종 출력은 공통 API v2의 14개 FITS 열을 사용합니다. 원본 기준 FITS의 두 열이나
이전 H/SZA 지도 형식과 구분하십시오.

원하는 좌표만 확인하는 예입니다. RA/DEC 입력은 degree, 날짜는 UTC입니다.
관측 조건 밖이면 밝기가 NaN이므로 `visible`도 확인합니다.

```python
from sdrift_sky import SkyModel

sky = SkyModel()
p = sky.point(
    30, 10, "2026-10-19T00:00:00Z", model="total", band="VIS",
    orbit="noon-midnight", orbital_phase_deg=90,
)
print(p["ab_mag_arcsec2"], p["intensity_mjy_sr"], p["visible"])
```

웹 화면 확인:

```bash
python3 -m http.server 8765 --directory docs
```

브라우저에서 `http://localhost:8765/`를 엽니다. 이미 켜 둔 Python/GUI는 새 모델로
재시작하고 브라우저도 새로고침하여 이전 데이터 캐시와 혼동하지 않도록 합니다.

## 10. 생성 결과 검증

### 10.1 원본·출력·환경 기록

생성 후에는 원본 archive, 코드 버전, Python 패키지 버전, 명령행, 로그와 함께
FITS/NPZ의 JSON metadata를 보관합니다.
SHA-256이 달라졌다고 물리 결과가 항상 달라진 것은 아니지만, 이전 데이터와
같은 입력이라고 간주해서도 안 됩니다.

```bash
python3 -m pip freeze > airglow_data/build_python_packages.txt
```

다음은 생성 파일의 존재, FITS 해상도·좌표계 및 exporter hash를 점검하는 예입니다.
비압축 배열의 hash를 사용한다는 점에 유의하십시오.

```python
from pathlib import Path
import gzip
import hashlib
import json
import numpy as np
from astropy.io import fits

root = Path(".")
for band in ["UV1", "UV2", "UV3", "VIS"]:
    for directory, prefix in [("zodi_data", "zodiacal"), ("dgl_data", "dgl")]:
        path = root / directory / f"{prefix}_{band}_nside128_template.fits.gz"
        with fits.open(path, checksum=True) as hd:
            hd.verify("exception")
            assert hd[1].header["NSIDE"] == 128
            assert hd[1].header["ORDERING"] == "RING"
            assert hd[1].header["COORDSYS"] == "C"
            values = hd[1].data["INTENSITY_MJY_SR"].ravel()
            assert values.size == 12 * 128**2
            assert np.all(np.isfinite(values) & (values >= 0))

manifest = json.loads((root / "docs/data/manifest.json").read_text())
assert set(manifest["bands"]) == {"UV1", "UV2", "UV3", "VIS"}
assert "dgl" in manifest
assert set(manifest["airglow"]) == {"glow", "combined"}

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

FITS에 CHECKSUM이 없으면 `checksum=True`만으로 체크섬 존재를 보장하지 않습니다.
위 코드는 FITS 구조와 값도 검사하며 변환된 배열은 SHA-256으로 검사합니다.

### 10.2 로컬 검증 스크립트

개발 폴더에 다음 테스트가 있으면 실행합니다.
Node.js를 쓰는 검증은 Python에서 JavaScript 구현도 호출합니다.

```bash
python3 -m unittest test_dgl test_airglow test_airglow_components
python3 test_web_sky.py
python3 -m unittest test_airglow_orbit test_sky_coordinates test_python_interface
```

LOS LUT를 만들었을 때만 추가로 실행합니다.

```bash
python3 -m unittest test_airglow_lut
```

이 검사는 수치 적분·단위·마스크·구현 간 일치를 확인하며, 모든 성분의 완전성이나
실제 위성 관측에 대한 물리 정확도를 보증하지 않습니다.

## 11. 입력을 바꿨을 때의 재생성 범위

| 바꾼 것 | 다시 만들 것 |
|---|---|
| Zodiacal 원본·기준 UTC·파장 표본 | Zodiacal FITS → 전체 export |
| IRIS/SFD·correlation·DGL 처방 | DGL FITS/metadata → DGL만 또는 전체 export |
| GLOW 원본·F10.7·Ap·SZA 간격·고도 상한 | GLOW NPZ → PALACE 결합 NPZ → 전체 export |
| PALACE 원본·층 폭·성분 선택·분자 분기 | 결합 NPZ → 전체 export |
| 위 Airglow 변경 후 이전 H/SZA용 LUT도 사용 | 해당 LUT와 지도 캐시도 재생성 |
| 관측 UTC·궤도·위상·H·ELV·회피각 | 원본 재생성 불필요; 최종 지도 재계산 |

공통 API가 읽는 `docs/data/`를 갱신하지 않으면 NPZ나 기준 FITS를 다시 만들어도
화면의 값은 바뀌지 않습니다. manifest와 대응 배열을 한 묶음으로 갱신합니다.

## 12. GitHub에 배포할 파일

이용자에게 필요한 것은 갱신한 `docs/data/`, 기존 웹 코드와 공통 Python 코드,
requirements 및 설명서입니다. 개발 폴더와 공개 clone을 나눈 현재 구성에서는
검증 후 `docs/data/`를 `web_publish/docs/data/`로 복사하여 배포합니다.
이 문서를 추가하는 작업만으로 데이터를 재생성하거나 다시 복사할 필요는 없습니다.

원본 GLOW/PALACE 소스, IRIS/SFD 대용량 FITS, 빌드 생성물, 이전 지도와 로그는
이용자 실행에 필요하지 않으므로 일반 사용자 배포에 포함하지 않습니다.
`.gitignore`는 기존 방침에 따라 push 대상에서 제외합니다.
재현용 전체 개발 환경은 공개 배포본과 별도로 백업합니다.

각 원본의 라이선스와 출처를 보존합니다. GLOW 계산의 이용은 GLOW academic
research license를 따르며, PALACE 코드는 GPLv3, 배포 모델 데이터는 CC BY 4.0입니다.
이 프로젝트는 GLOW 등을 이용한 계산 결과를 표시하는 소프트웨어입니다.

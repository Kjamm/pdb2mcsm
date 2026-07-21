# pdb2mcsm

[mCSM-PPI2](https://biosig.lab.uq.edu.au/mcsm_ppi2/) ΔΔG 분석을 위해 단백질 복합체 PDB 파일을 정리해주는 도구입니다.

MD 시뮬레이션 출력이나 그 밖의 raw 구조 파일은 mCSM-PPI2에 그대로 못 넣는 경우가 많습니다. 사슬이 합쳐져 있거나 chain 표시가 없고, 잔기 번호가 1부터 다시 매겨져 있고, 물·이온이 가득 붙어 있고, 히스티딘이 `HIE` 같은 MD 전용 이름으로 적혀 있기 때문입니다. `pdb2mcsm`은 이런 문제를 자동으로 고쳐서 서버가 받아들이는 형태로 만들어줍니다.

## 코랩에서 바로 쓰기 (설치 불필요)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1is8nIMIMRTz3BbJSpFM4oRquwRRakQ8I?usp=sharing)

PDB를 업로드하면 정리된 `_ready.pdb`를 받고, 원하면 alanine-scanning 돌연변이 목록까지 생성할 수 있습니다. 전부 입력칸으로 조작하며 코딩이 필요 없습니다.

## 기능

**자동 정리**
- **사슬 분리** — 우선순위: 기존 chain ID → `TER` 레코드 → 잔기 번호 점프 (추측한 경우 경고 표시).
- **사슬 개수 무관** — 2개, 3개, 항체 복합체 등 모두 처리.
- **물·이온·결정화 첨가물 제거** (HOH, NA, CL, SO4, GOL 등).
- **수소 제거** (mCSM이 다시 붙임).
- **MD 잔기명 정규화** — HIE/HID/HIP → HIS, CYX → CYS, MSE → MET 등.
- **인산화 잔기 처리** — SEP/S1P → SER, TPO/T1P → THR, PTR/Y1P → TYR (인산기 원자 제거 포함).
- **삽입 코드 보존** (예: 100A).
- 깨끗한 chain ID 부여 (원본이 유효하면 유지).

**선택 기능**
- 사슬 **번호 이동** — 고정 오프셋으로 (예: 파일 번호를 논문 번호로 맞추기).
- **alanine-scanning 돌연변이 목록** 생성 (mCSM-PPI2 배치 형식).

## 빠른 시작 (Python)

```python
import core

# 1. 읽고 사슬 분리
header, chains, warnings, method = core.parse_pdb("내복합체.pdb")
core.assign_chain_ids(chains)

# 2. 검토 (꼭 확인!)
print("분리 방법:", method)
for row in core.summarize(chains):
    print(row)          # (체인, 잔기수, 첫잔기, 끝잔기, 분리근거)
for w in warnings:
    print("경고:", w)

# 3. (선택) 사슬 번호 이동: 파일 36 → 논문 98 이면 offset 62
core.renumber_chain(chains, "B", 62)

# 4. 정리된 파일 저장 → 이 파일을 mCSM-PPI2에 업로드
core.write_clean("내복합체_ready.pdb", header, chains)

# 5. (선택) 배치용 alanine-scanning 목록
print("\n".join(core.build_ala_scan(chains, "B", 95, 110)))
```

## 이후 mCSM-PPI2에서

1. mCSM-PPI2의 배치(다중 돌연변이) 페이지를 엽니다.
2. `내복합체_ready.pdb`를 업로드합니다.
3. 돌연변이시킬 사슬을 선택합니다.
4. alanine-scanning 목록을 붙여넣습니다.
5. 실행 → 잔기별 ΔΔG. 크게 음수인 잔기가 결합에 중요한 hotspot입니다.

## 함수 목록

| 함수 | 역할 |
|---|---|
| `parse_pdb(path)` | PDB를 읽어 사슬 분리; `(header, chains, warnings, method)` 반환 |
| `assign_chain_ids(chains)` | 각 사슬에 깨끗한 A/B/C ID 부여 |
| `write_clean(path, header, chains)` | 정리된 PDB 저장 (수소·인산기 원자 제거) |
| `renumber_chain(chains, chain_id, offset)` | 사슬 잔기 번호를 `offset`만큼 이동 |
| `build_ala_scan(chains, chain_id, start, end)` | mCSM-PPI2 alanine-scanning 목록 생성 |
| `summarize(chains)` | 검토용 사슬별 요약 |
| `read_atoms(path)` | 저수준: ATOM/HETATM을 `Atom` 객체로 파싱 |

## 한계

`pdb2mcsm`은 파일의 **구조·정리 문제**를 해결합니다. 다음은 **처리하지 않습니다**:

- **자동 번호 매핑** — 오프셋을 직접 알아야 합니다 (`renumber_chain`은 고정 이동).
- **인산화 외 다른 변형 잔기** — 아세틸화, 메틸화, 당화 등은 미포함. 필요하면 `NONSTD_MAP`에 추가하세요.
- **리간드·보조인자** (약물, 헴, ATP, 금속 클러스터) — 물·이온만 제거합니다.
- **구조 결함** — 빠진 원자, 끊긴 백본, 원자 충돌, 대체 위치(altloc)는 고치지 않습니다.
- **핵산** — 단백질-DNA/RNA 복합체는 대상 밖 (mCSM-PPI2는 단백질-단백질 전용).

일부 한계는 이 도구가 아니라 **mCSM-PPI2 자체의 제약**입니다: 잔기를 표준형으로 바꾸는 순간 인산화 *효과*는 사라집니다. 단백질-리간드·단백질-핵산 결합은 지원되지 않으며, 매우 큰 복합체는 서버 크기 제한을 넘을 수 있습니다.

## 참고

- 기본값은 **원본 잔기 번호 유지**입니다. mCSM에 넣는 번호는 정리된 파일에 적힌 번호와 같아야 합니다. `summarize()`나 `chains[i].residues`로 확인하세요.
- 인산화 잔기를 표준형으로 바꾸면 모델에서 인산화가 사라집니다. 인산화 자체가 결합에 중요한 연구라면 결과 해석에 주의하거나, 변형 잔기를 지원하는 방법을 쓰세요.

## 라이선스

MIT

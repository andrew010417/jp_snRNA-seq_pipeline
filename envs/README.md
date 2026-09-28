# 실행 환경

## 필요한 것

| | |
|---|---|
| R | 4.3 이상 (검증: **4.4.3**) |
| 메모리 | 세포 10만 개 기준 약 32 GB |

## 설치

```bash
Rscript envs/install.R
```

## 패키지

| 패키지 | 쓰는 곳 | 필수 |
|---|---|---|
| Seurat (5.x) | 전 단계 | ✅ |
| DoubletFinder | 02 더블렛 | ✅ |
| scran · scuttle · batchelor | 04 정규화(scran 방식) | ✅ |
| data.table · yaml · ggplot2 · patchwork | 공통 | ✅ |
| R.utils | 가짜 데이터 생성 | 테스트용 |
| presto | 08 마커 (약 20배 빠름) | 선택 |
| harmony | 06 배치 보정 | **선택** |
| clustree | 07 해상도 고르기 | 선택 |

### harmony 가 없다면

`config/config.yaml` 에서 `batch.method` 를 **`rpca`** 로 바꾸세요.
Seurat 내장 기능이라 추가 설치 없이 돕니다.

## 알려진 버전 충돌

| 증상 | 원인 | 영향 |
|---|---|---|
| `Unknown guide: edge_colourbar` | `clustree` 가 **ggplot2 4.x** 와 맞지 않음 | **clustree 그림 한 장만** 안 나오고 파이프라인은 계속 진행됩니다 |

그림 저장은 전부 `tryCatch` 로 감싸 두었습니다.
**그림이 실패해도 계산 결과는 잃지 않습니다.** 실패 사유는 `run/logs/plot_failed_*.txt` 에 남습니다.

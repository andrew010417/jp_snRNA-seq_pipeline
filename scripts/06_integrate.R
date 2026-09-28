#!/usr/bin/env Rscript
## 06 — 샘플·환자·실험일 차이(배치 효과)를 지웁니다.
##   ⚠️ 이 단계는 PCA 뒤, 이웃 그래프 앞이어야 합니다.
##      보정 전 좌표로 이웃을 정하면 배치 효과가 클러스터에 그대로 남습니다.
##   harmony : PCA 좌표를 배치별로 밀어 맞춤 (가장 흔함)
##   rpca    : Seurat 내장. 추가 패키지 없이 됨
##   none    : 배치가 하나거나 보정이 필요 없을 때
## 출력: run/rds/integrated.rds  (reduction 이름은 항상 "integrated")
suppressPackageStartupMessages({ library(Seurat) })
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_common.R"))

obj    <- readRDS(file.path(RDS, "pca.rds"))
method <- CFG$batch$method
key    <- CFG$batch$key
nb     <- length(unique(obj[[key]][, 1]))
say("배치 보정: ", method, "  (", key, " ", nb, "개)")

if (nb < 2 || identical(method, "none")) {
  say("보정하지 않고 PCA 좌표를 그대로 씁니다")
  obj[["integrated"]] <- obj[["pca"]]
} else if (identical(method, "harmony")) {
  if (!requireNamespace("harmony", quietly = TRUE))
    stop("harmony 패키지가 없습니다. envs/install.R 을 실행하거나 batch.method 를 'rpca'로 바꾸세요.")
  obj <- harmony::RunHarmony(obj, group.by.vars = key,
                             reduction.use = "pca", reduction.save = "integrated",
                             verbose = FALSE)
} else if (identical(method, "rpca")) {
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj[[key]][, 1])
  obj <- IntegrateLayers(obj, method = RPCAIntegration,
                         orig.reduction = "pca", new.reduction = "integrated",
                         verbose = FALSE)
  obj <- JoinLayers(obj)
} else stop("batch.method 는 harmony | rpca | none 중 하나여야 합니다: ", method)

saveRDS(obj, file.path(RDS, "integrated.rds"))
say("저장: integrated.rds")
say("완료: 06_integrate")

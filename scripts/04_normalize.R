#!/usr/bin/env Rscript
## 04 — 핵마다 읽힌 총량이 다르므로 "같은 깊이로 읽었다면"의 값으로 맞춥니다.
##   scran        : 비슷한 핵끼리 묶어 size factor를 추정 → 희소한 데이터에 강함
##   lognormalize : Seurat 기본. 빠름
## 출력: run/rds/normalized.rds
suppressPackageStartupMessages({ library(Seurat) })
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_common.R"))

obj <- readRDS(file.path(RDS, "merged_filtered.rds"))
method <- CFG$normalization$method
say("정규화 방법: ", method)

if (identical(method, "scran")) {
  suppressPackageStartupMessages({
    library(scran); library(scuttle); library(SingleCellExperiment); library(batchelor)
  })
  sce <- SingleCellExperiment(list(counts = GetAssayData(obj, layer = "counts")))
  sce$sample <- obj$sample

  say("quickCluster…")
  cl <- scran::quickCluster(sce, min.size = min(100, floor(ncol(sce) / 4)))
  say("computeSumFactors…")
  sce <- scran::computeSumFactors(sce, clusters = cl, min.mean = 0.1)

  ## 샘플이 여러 개면 배치 간 규모까지 맞춥니다
  if (length(unique(sce$sample)) > 1) {
    say("multiBatchNorm…")
    sce <- batchelor::multiBatchNorm(sce, batch = sce$sample)
  } else {
    sce <- scuttle::logNormCounts(sce)
  }
  LayerData(obj, layer = "data") <- as(logcounts(sce), "dgCMatrix")
  obj$size_factor <- sizeFactors(sce)
} else {
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 1e4, verbose = FALSE)
}

saveRDS(obj, file.path(RDS, "normalized.rds"))
say("저장: normalized.rds")
say("완료: 04_normalize")

#!/usr/bin/env Rscript
## 05 — 변이가 큰 유전자만 고르고(HVG), 분산을 맞춘 뒤(Scaling), 차원을 줄입니다(PCA).
## 출력: run/rds/pca.rds, run/figures/elbow.png, run/tables/hvg.txt
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_common.R"))

obj <- readRDS(file.path(RDS, "normalized.rds"))

n_hvg <- min(CFG$features$n_hvg, nrow(obj))
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = n_hvg, verbose = FALSE)
writeLines(VariableFeatures(obj), file.path(TAB, "hvg.txt"))
say("HVG ", length(VariableFeatures(obj)), "개 선택")

vr <- CFG$features$vars_to_regress
obj <- ScaleData(obj, vars.to.regress = if (length(vr)) unlist(vr) else NULL, verbose = FALSE)

npcs <- min(CFG$reduction$n_pcs, ncol(obj) - 1, length(VariableFeatures(obj)) - 1)
obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
say("PCA ", npcs, "차원")

## Elbow plot — 몇 개의 PC를 쓸지 눈으로 정하기 위한 그림
save_plot(ElbowPlot(obj, ndims = npcs) + theme_minimal(base_size = 11), "elbow", w = 7, h = 5)

saveRDS(obj, file.path(RDS, "pca.rds"))
say("저장: pca.rds")
say("완료: 05_hvg_scale_pca")

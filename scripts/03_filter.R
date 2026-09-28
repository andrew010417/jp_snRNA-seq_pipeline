#!/usr/bin/env Rscript
## 03 — QC 기준을 적용해 핵을 거르고, 거의 안 잡히는 유전자를 뺀 뒤 하나로 합칩니다.
## 출력: run/rds/merged_filtered.rds, run/qc/qc_cascade.tsv
suppressPackageStartupMessages({ library(Seurat) })
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_common.R"))

samples <- sample_list()
q <- CFG$qc
cascade <- list(); objs <- list()

for (s in samples) {
  f_df <- file.path(RDS, paste0(s, "_df.rds"))
  o <- readRDS(if (file.exists(f_df)) f_df else file.path(RDS, paste0(s, "_raw.rds")))
  n0 <- ncol(o)

  ## 걸러낸 뒤 남은 핵이 하나도 없으면 원인을 짚어 멈춥니다
  guard <- function(keep, what, hint) {
    if (sum(keep) == 0)
      stop(sprintf("[%s] %s 단계에서 모든 핵이 걸러졌습니다.\n  %s\n  run/qc/qc_metrics.tsv 의 중앙값을 보고 config/config.yaml 의 qc 항목을 조정하세요.",
                   s, what, hint), call. = FALSE)
    keep
  }

  ## 1) 더블렛
  if (!is.null(o$doublet)) o <- o[, guard(!o$doublet, "더블렛 제거", "모든 핵이 더블렛으로 판정됐습니다. doublet.fixed_rate 를 확인하세요.")]
  n1 <- ncol(o)

  ## 2) 핵 QC — 고정 임계값 또는 샘플별 중앙값 비율
  if (isTRUE(q$use_median_ratio)) {
    keep <- o$nFeature_RNA > q$median_ratio * median(o$nFeature_RNA) &
            o$nCount_RNA   > q$median_ratio * median(o$nCount_RNA)
  } else {
    keep <- o$nFeature_RNA >= q$min_features & o$nCount_RNA >= q$min_counts
  }
  o <- o[, guard(keep, "핵 QC",
                 sprintf("기준: min_features=%s, min_counts=%s / 이 샘플 중앙값: nFeature=%s, nCount=%s",
                         q$min_features, q$min_counts,
                         median(o$nFeature_RNA), median(o$nCount_RNA)))]
  n2 <- ncol(o)

  ## 3) 미토콘드리아 — snRNA에서는 세포질 오염 지표입니다
  o <- o[, guard(o$percent_mt < q$max_percent_mt, "미토콘드리아 필터",
                 sprintf("기준: percent_mt < %s / 이 샘플 중앙값: %.2f%%",
                         q$max_percent_mt, median(o$percent_mt)))]
  n3 <- ncol(o)

  cascade[[s]] <- data.table(sample = s, 시작 = n0, 더블렛제거후 = n1,
                             핵QC후 = n2, 미토필터후 = n3,
                             최종유지율 = round(n3 / n0 * 100, 1))
  objs[[s]] <- o
  say(s, ": ", n0, " → ", n3, " (", round(n3 / n0 * 100, 1), "%)")
}

cascade <- rbindlist(cascade)
fwrite(cascade, file.path(QCD, "qc_cascade.tsv"), sep = "\t")
print(cascade)

## 합치기
merged <- if (length(objs) == 1) objs[[1]] else
  merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))
merged <- JoinLayers(merged)

## 유전자 필터 — 너무 적은 핵에서만 검출된 유전자 제거
cnt <- GetAssayData(merged, layer = "counts")
ncell_per_gene <- Matrix::rowSums(cnt > 0)
keep_g <- ncell_per_gene >= CFG$qc$min_cells_per_gene
if (sum(keep_g) == 0)
  stop("min_cells_per_gene 기준으로 모든 유전자가 제거됐습니다. 값을 낮추세요.", call. = FALSE)
say("유전자 ", nrow(merged), " → ", sum(keep_g))
merged <- subset(merged, features = rownames(merged)[keep_g])

saveRDS(merged, file.path(RDS, "merged_filtered.rds"))
say("저장: merged_filtered.rds  (핵 ", ncol(merged), " · 유전자 ", nrow(merged), ")")
say("완료: 03_filter")

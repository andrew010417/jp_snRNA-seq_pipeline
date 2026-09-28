#!/usr/bin/env Rscript
## 01 — 샘플별 10x 행렬을 읽어 Seurat 객체로 만들고 QC 지표를 계산합니다.
## 출력: run/rds/<sample>_raw.rds, run/qc/qc_metrics.tsv, run/figures/qc_violin.png
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_common.R"))

samples <- sample_list()
say("샘플 ", length(samples), "개: ", paste(samples, collapse = ", "))

metrics <- rbindlist(lapply(samples, function(s) {
  out <- file.path(RDS, paste0(s, "_raw.rds"))
  obj <- checkpoint(out, {
    say(s, " 읽는 중…")
    mtx <- Read10X(data.dir = file.path(IN, s))
    o <- CreateSeuratObject(counts = mtx, project = s,
                            min.cells = 0, min.features = 0)
    o$sample <- s
    o$percent_mt <- PercentageFeatureSet(o, pattern = MT_PATTERN)
    o
  })
  data.table(sample = s, n_cells = ncol(obj), n_genes = nrow(obj),
             median_counts = median(obj$nCount_RNA),
             median_features = median(obj$nFeature_RNA),
             median_percent_mt = median(obj$percent_mt),
             pct_mt_over_5 = mean(obj$percent_mt > 5) * 100)
}))

fwrite(metrics, file.path(QCD, "qc_metrics.tsv"), sep = "\t")
say("QC 지표 저장: qc/qc_metrics.tsv")
print(metrics)

## QC 분포 그림 — 임계값을 눈으로 정하기 위한 것
md <- rbindlist(lapply(samples, function(s) {
  o <- readRDS(file.path(RDS, paste0(s, "_raw.rds")))
  data.table(sample = s, nCount_RNA = o$nCount_RNA,
             nFeature_RNA = o$nFeature_RNA, percent_mt = o$percent_mt)
}))
long <- melt(md, id.vars = "sample")
p <- ggplot(long, aes(sample, value)) +
  geom_violin(scale = "width", fill = "#4C72B0", alpha = .6, colour = NA) +
  facet_wrap(~ variable, scales = "free_y", ncol = 1) +
  scale_y_continuous(trans = "log1p") +
  labs(x = NULL, y = NULL, title = "QC 분포 (y축 log1p)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, "qc_violin", w = max(6, length(samples) * 0.8), h = 8)
say("완료: 01_load_qc")

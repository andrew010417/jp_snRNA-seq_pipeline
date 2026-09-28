#!/usr/bin/env Rscript
## 08 — 클러스터마다 특징 유전자를 뽑습니다. 이걸 보고 세포 유형 이름을 붙입니다.
## 출력: run/tables/markers_all.tsv, markers_top.tsv, run/figures/marker_dotplot.png
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_common.R"))

obj <- readRDS(file.path(RDS, "clustered.rds"))
say("클러스터 ", nlevels(factor(Idents(obj))), "개에서 마커 찾는 중…")

mk <- FindAllMarkers(obj, only.pos = TRUE,
                     logfc.threshold = CFG$markers$logfc_threshold,
                     min.pct = CFG$markers$min_pct, verbose = FALSE)
if (!nrow(mk)) { say("마커가 없습니다 — 임계값을 낮춰 보세요"); quit(save = "no") }

mk <- as.data.table(mk)
fwrite(mk, file.path(TAB, "markers_all.tsv"), sep = "\t")

top <- mk[order(cluster, -avg_log2FC)][, head(.SD, CFG$markers$top_n), by = cluster]
fwrite(top, file.path(TAB, "markers_top.tsv"), sep = "\t")
say("마커 저장: tables/markers_all.tsv · markers_top.tsv")

top3 <- unique(mk[order(cluster, -avg_log2FC)][, head(.SD, 3), by = cluster]$gene)
if (length(top3) > 1) {
  p <- DotPlot(obj, features = top3) + RotatedAxis() +
    theme_minimal(base_size = 10) + labs(title = "클러스터별 상위 마커")
  save_plot(p, "marker_dotplot", w = max(8, length(top3) * 0.35), h = 6)
}

say("")
say("다음 단계: tables/markers_top.tsv 를 보고 클러스터마다 세포 유형 이름을 붙이세요.")
say("  자동 주석 도구를 쓰려면 SingleR · CellTypist · Azimuth 를 검토하세요.")
say("완료: 08_markers")

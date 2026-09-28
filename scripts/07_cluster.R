#!/usr/bin/env Rscript
## 07 — 이웃 그래프를 만들고, 비슷한 핵끼리 묶고(클러스터링), 2차원으로 펼칩니다(UMAP).
##   ⚠️ UMAP은 클러스터링 결과물이 아니라 점검용 그림입니다.
##      UMAP 위의 거리를 해석하지 마세요.
## 출력: run/rds/clustered.rds, run/figures/umap_*.png, run/figures/clustree.png
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_common.R"))

obj  <- readRDS(file.path(RDS, "integrated.rds"))
dims <- seq_len(min(CFG$reduction$n_pcs, ncol(Embeddings(obj, "integrated"))))

obj <- FindNeighbors(obj, reduction = "integrated", dims = dims,
                     k.param = CFG$cluster$k_param, verbose = FALSE)

alg  <- if (identical(CFG$cluster$algorithm, "leiden")) 4L else 1L
res  <- unlist(CFG$cluster$resolutions)
obj  <- FindClusters(obj, resolution = res, algorithm = alg,
                     random.seed = CFG$project$seed, verbose = FALSE)

## 해상도를 올릴 때 세포가 어디서 갈라지는지 — 해상도 고르기용
if (requireNamespace("clustree", quietly = TRUE) && length(res) > 1) {
  p <- clustree::clustree(obj, prefix = "RNA_snn_res.")
  save_plot(p, "clustree", w = 9, h = 8)
}

sel <- paste0("RNA_snn_res.", CFG$cluster$select)
if (!sel %in% colnames(obj@meta.data))
  stop("cluster.select 값이 resolutions 안에 없습니다: ", CFG$cluster$select)
Idents(obj) <- obj$seurat_clusters <- obj@meta.data[[sel]]
say("해상도 ", CFG$cluster$select, " → 클러스터 ", nlevels(factor(Idents(obj))), "개")

obj <- RunUMAP(obj, reduction = "integrated", dims = dims,
               n.neighbors = CFG$umap$n_neighbors, min.dist = CFG$umap$min_dist,
               metric = CFG$umap$metric, seed.use = CFG$project$seed, verbose = FALSE)

save_plot(DimPlot(obj, label = TRUE, repel = TRUE) + theme_minimal(base_size = 11) +
            labs(title = paste0("클러스터 (resolution ", CFG$cluster$select, ")")),
          "umap_clusters")
save_plot(DimPlot(obj, group.by = CFG$batch$key) + theme_minimal(base_size = 11) +
            labs(title = "배치별 — 잘 섞였는지 점검"),
          "umap_batch")

fwrite(as.data.table(table(sample = obj$sample, cluster = Idents(obj))),
       file.path(TAB, "cluster_by_sample.tsv"), sep = "\t")

saveRDS(obj, file.path(RDS, "clustered.rds"))
say("저장: clustered.rds")
say("완료: 07_cluster")

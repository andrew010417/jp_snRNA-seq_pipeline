#!/usr/bin/env Rscript
## 02 — DoubletFinder로 샘플별 더블렛(한 방울에 핵 2개)을 찾습니다.
## 샘플 하나씩 도는 것이 안전합니다:  Rscript 02_doublet.R <sample>
## 인자가 없으면 모든 샘플을 차례로 돕니다.
## 출력: run/rds/<sample>_df.rds, run/qc/doublet_rates.tsv
suppressPackageStartupMessages({ library(Seurat); library(DoubletFinder) })
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), "00_common.R"))

if (!isTRUE(CFG$doublet$enabled)) { say("doublet.enabled = false — 건너뜁니다"); quit(save = "no") }

args <- commandArgs(TRUE)
targets <- if (length(args)) args else sample_list()

rates <- rbindlist(lapply(targets, function(s) {
  out <- file.path(RDS, paste0(s, "_df.rds"))
  if (file.exists(out)) { say("건너뜀: ", s); o <- readRDS(out) }
  else {
    say(s, " DoubletFinder 실행 중…")
    o <- readRDS(file.path(RDS, paste0(s, "_raw.rds")))
    N <- ncol(o)
    pcs <- seq_len(min(CFG$doublet$pcs, N - 2))

    o <- NormalizeData(o, verbose = FALSE)
    o <- FindVariableFeatures(o, selection.method = "vst",
                              nfeatures = min(CFG$features$n_hvg, nrow(o)), verbose = FALSE)
    o <- ScaleData(o, verbose = FALSE)
    o <- RunPCA(o, npcs = min(30, N - 1), verbose = FALSE)

    ## pK(이웃 크기) 탐색 — 데이터마다 최적값이 달라 반드시 sweep 합니다
    sw   <- paramSweep(o, PCs = pcs, sct = FALSE, num.cores = 1)
    bcm  <- find.pK(summarizeSweep(sw, GT = FALSE))
    bcm$pK <- as.numeric(as.character(bcm$pK))
    pK_opt <- bcm$pK[which.max(bcm$BCmetric)]

    rate <- if (is.null(CFG$doublet$fixed_rate)) tenx_rate(N) else CFG$doublet$fixed_rate
    nExp <- round(rate * N)

    o <- doubletFinder(o, PCs = pcs, pN = CFG$doublet$pN, pK = pK_opt,
                       nExp = nExp, sct = FALSE)
    cls <- grep("^DF.classifications", colnames(o@meta.data), value = TRUE)[1]
    o$doublet <- o@meta.data[[cls]] == "Doublet"
    o@meta.data$pK_used <- pK_opt
    saveRDS(o, out); say("저장: ", basename(out))
  }
  data.table(sample = s, n_cells = ncol(o), n_doublet = sum(o$doublet),
             pct_doublet = round(mean(o$doublet) * 100, 2),
             pK = o$pK_used[1])
}))

f <- file.path(QCD, "doublet_rates.tsv")
fwrite(rates, f, sep = "\t")
say("더블렛 비율 저장: qc/doublet_rates.tsv")
print(rates)
say("완료: 02_doublet")

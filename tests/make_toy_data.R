#!/usr/bin/env Rscript
## 파이프라인이 끝까지 도는지 확인하기 위한 작은 가짜 데이터를 만듭니다.
## 실제 생물학적 의미는 없습니다. 배선 점검용입니다.
##   사용법: Rscript tests/make_toy_data.R <출력폴더> [샘플수] [샘플당세포수] [유전자수]
suppressPackageStartupMessages({ library(Matrix) })
a <- commandArgs(TRUE)
out      <- if (length(a) >= 1) a[1] else "data/raw"
n_sample <- if (length(a) >= 2) as.integer(a[2]) else 3
n_cell   <- if (length(a) >= 3) as.integer(a[3]) else 300
n_gene   <- if (length(a) >= 4) as.integer(a[4]) else 400
set.seed(1)

## 유전자 이름 — 미토콘드리아 10개를 섞어 percent_mt 계산이 되는지 확인
genes <- c(sprintf("MT-G%02d", 1:10), sprintf("GENE%04d", seq_len(n_gene - 10)))
n_type <- 4                                    # 가상의 세포 유형 4종
type_mu <- matrix(rgamma(n_type * n_gene, shape = .4, scale = 2),
                  nrow = n_type, dimnames = list(NULL, genes))
for (k in seq_len(n_type)) {                   # 유형마다 마커 30개를 크게 올림
  idx <- sample(11:n_gene, 30); type_mu[k, idx] <- type_mu[k, idx] + 12
}

for (s in seq_len(n_sample)) {
  nm  <- sprintf("SAMPLE_%s", LETTERS[s])
  dir.create(file.path(out, nm), recursive = TRUE, showWarnings = FALSE)
  tp     <- sample(seq_len(n_type), n_cell, replace = TRUE)
  depth  <- runif(n_cell, 0.7, 1.4)            # 핵마다 다른 시퀀싱 깊이
  batch  <- rgamma(n_gene, shape = 8, rate = 8)  # 샘플 고유 배치 효과
  m <- vapply(seq_len(n_cell), function(i)
        rpois(n_gene, type_mu[tp[i], ] * depth[i] * batch), numeric(n_gene))
  rownames(m) <- genes
  colnames(m) <- sprintf("%s-1", replicate(n_cell,
                    paste(sample(c("A","C","G","T"), 16, TRUE), collapse = "")))
  m <- as(as(m, "dgCMatrix"), "generalMatrix")

  writeMM(m, file.path(out, nm, "matrix.mtx"))
  R.utils::gzip(file.path(out, nm, "matrix.mtx"), overwrite = TRUE)
  writeLines(colnames(m), gzfile(file.path(out, nm, "barcodes.tsv.gz")))
  writeLines(paste(genes, genes, "Gene Expression", sep = "\t"),
             gzfile(file.path(out, nm, "features.tsv.gz")))
  cat("만듦:", nm, "-", n_gene, "유전자 x", n_cell, "세포\n")
}
cat("완료 →", out, "\n")

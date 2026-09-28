#!/usr/bin/env Rscript
## 필요한 R 패키지를 설치합니다.
##   Rscript envs/install.R
## 이미 있는 패키지는 건너뜁니다.

cran <- c("Matrix", "data.table", "yaml", "ggplot2", "patchwork",
          "Seurat", "remotes", "R.utils")
bioc <- c("scran", "scuttle", "batchelor", "SingleCellExperiment")
opt  <- c("harmony", "clustree")            # 없어도 파이프라인은 돕니다
gh   <- c(DoubletFinder = "chris-mcginnis-ucsf/DoubletFinder",
          presto        = "immunogenomics/presto")

need <- function(p) !requireNamespace(p, quietly = TRUE)

miss <- Filter(need, c(cran, opt))
if (length(miss)) { cat("CRAN 설치:", paste(miss, collapse = ", "), "\n")
                    install.packages(miss, repos = "https://cloud.r-project.org") }

if (any(sapply(bioc, need))) {
  if (need("BiocManager")) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  m <- Filter(need, bioc)
  cat("Bioconductor 설치:", paste(m, collapse = ", "), "\n")
  BiocManager::install(m, ask = FALSE, update = FALSE)
}

for (p in names(gh)) if (need(p)) {
  cat("GitHub 설치:", p, "\n"); remotes::install_github(gh[[p]], upgrade = "never")
}

cat("\n설치 상태\n")
for (p in c(cran, bioc, opt, names(gh)))
  cat(sprintf("  %-22s %s\n", p,
      if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "없음"))
cat("\n'없음'이 harmony·clustree 뿐이면 그대로 쓰셔도 됩니다",
    "(batch.method 를 'rpca' 로 두면 harmony 없이 돕니다).\n")

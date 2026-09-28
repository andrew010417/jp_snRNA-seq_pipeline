## 00_common.R — 모든 단계가 공유하는 설정·경로·도우미
suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(yaml)
})

## ---- 설정 읽기 -------------------------------------------------------------
## 환경변수 SNRNA_CONFIG 로 다른 설정 파일을 지정할 수 있습니다.
find_root <- function() {
  a <- grep("--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(normalizePath(file.path(dirname(sub("--file=", "", a[1])), "..")))
  normalizePath("..")
}
ROOT <- Sys.getenv("SNRNA_ROOT", unset = find_root())
CFG_PATH <- Sys.getenv("SNRNA_CONFIG", unset = file.path(ROOT, "config/config.yaml"))
stopifnot(file.exists(CFG_PATH))
CFG <- yaml::read_yaml(CFG_PATH)

abspath <- function(p) if (startsWith(p, "/")) p else file.path(ROOT, p)
IN   <- abspath(CFG$paths$input_dir)
OUT  <- abspath(CFG$paths$out_dir)
RDS  <- file.path(OUT, "rds")
QCD  <- file.path(OUT, "qc")
FIG  <- file.path(OUT, "figures")
TAB  <- file.path(OUT, "tables")
LOG  <- file.path(OUT, "logs")
for (d in c(RDS, QCD, FIG, TAB, LOG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

set.seed(CFG$project$seed)

MT_PATTERN <- if (identical(CFG$project$species, "mouse")) "^mt-" else "^MT-"

## ---- 로그 ------------------------------------------------------------------
say <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## ---- 샘플 목록 -------------------------------------------------------------
## input_dir 아래의 각 폴더를 샘플 하나로 봅니다.
sample_list <- function() {
  d <- list.dirs(IN, recursive = FALSE, full.names = FALSE)
  d <- d[nzchar(d)]
  if (!length(d)) stop("입력 폴더에 샘플 폴더가 없습니다: ", IN)
  sort(d)
}

## ---- 10x Chromium 3' v3 기대 더블렛 비율 표 --------------------------------
## 출처: 10x Genomics Chromium Next GEM Single Cell 3' v3.1 User Guide
TENX_RATE <- data.frame(
  recovered = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000),
  rate      = c(0.004, 0.008, 0.016, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076))

tenx_rate <- function(n) {
  tb <- TENX_RATE
  vapply(n, function(x) {
    if (x <= tb$recovered[1]) return(tb$rate[1] * x / tb$recovered[1])
    if (x >= tb$recovered[nrow(tb)]) {
      sl <- diff(tail(tb$rate, 2)) / diff(tail(tb$recovered, 2))
      return(tb$rate[nrow(tb)] + sl * (x - tb$recovered[nrow(tb)]))
    }
    approx(tb$recovered, tb$rate, xout = x)$y
  }, numeric(1))
}

## ---- 체크포인트 ------------------------------------------------------------
## 이미 만들어진 산출물이 있으면 건너뜁니다. 지우고 다시 돌리면 재계산합니다.
checkpoint <- function(path, expr) {
  if (file.exists(path)) { say("건너뜀 (이미 있음): ", basename(path)); return(readRDS(path)) }
  v <- force(expr); saveRDS(v, path); say("저장: ", basename(path)); v
}

## ---- 그림 저장 -------------------------------------------------------------
## 그림은 보조 산출물입니다. 패키지 버전 충돌로 한 장이 실패해도
## 계산 결과를 잃지 않도록, 경고만 남기고 다음으로 넘어갑니다.
save_plot <- function(p, name, w = 9, h = 6) {
  f <- file.path(FIG, paste0(name, ".png"))
  tryCatch({
    ggplot2::ggsave(f, p, width = w, height = h, dpi = 200)
    say("그림 저장: figures/", basename(f))
  }, error = function(e) {
    ## 실패하면 반쯤 쓰다 만 파일이 남을 수 있으므로 지웁니다
    if (file.exists(f)) unlink(f)
    say("[경고] 그림 '", name, "' 을(를) 그리지 못했습니다 — 계속 진행합니다")
    say("        사유: ", conditionMessage(e))
    writeLines(conditionMessage(e), file.path(LOG, paste0("plot_failed_", name, ".txt")))
  })
  invisible(p)
}

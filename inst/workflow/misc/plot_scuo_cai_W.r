### This script plot X-Y protein production rates.

rm(list = ls())

suppressMessages(library(cubfits, quietly = TRUE))

source("00-set_env.r")

# pre processed phi.Obs
fn.in <- paste(prefix$data, "pre_process.rda", sep = "")
load(fn.in)

# w <- seqinr::caitab$sc
 w <- seqinr::caitab$ec
names(w) <- codon.low2up(rownames(seqinr::caitab))
CAI.W <- calc_cai_values(y, y.list, w = w)

# plot
# x-axis: predicted, y-axis: observed.
fn.in <- file.data$tsv
tmp <- log10(phi.Obs / mean(phi.Obs))

fn.out <- paste(prefix$plot, "prxy_scuo_cai_W.pdf", sep = "")
pdf(fn.out, width = 12, height = 4)
par(mfrow = c(1, 3))

plotprxy(CAI$cai, tmp, log10.x = FALSE, log10.y = FALSE,
         xlab = "CAI", ylab = "Observed Production Rate (log10)",
         main = "CAI vs phi.Obs")

plotprxy(CAI.W$cai, tmp, log10.x = FALSE, log10.y = FALSE,
         xlab = "CAI.W", ylab = "Observed Production Rate (log10)",
         main = "CAI.W vs phi.Obs")

plotprxy(CAI.W$cai, CAI$cai, log10.x = FALSE, log10.y = FALSE,
         xlab = "CAI.W", ylab = "CAI",
         main = "CAI.W vs CAI")

dev.off()


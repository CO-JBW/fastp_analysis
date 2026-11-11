# Draws QC figures and writes a metrics table from a fastp JSON report.
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
rep <- fromJSON(args[1])
sample <- args[2]
figdir <- "/results/figures"
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

open_png <- function(name, width = 1000) {
  png(file.path(figdir, name), width = width, height = 600, res = 110)
}
reads <- c("read1", "read2")

# mean quality per cycle, before and after filtering
open_png("fig1_quality_per_cycle.png", 1600)
par(mfrow = c(1, 2))
for (r in reads) {
  before <- rep[[paste0(r, "_before_filtering")]]$quality_curves$mean
  after <- rep[[paste0(r, "_after_filtering")]]$quality_curves$mean
  plot(before, type = "l", col = "red", lwd = 2, ylim = c(0, 42),
       xlab = "cycle (bp)", ylab = "mean Phred score",
       main = paste(r, "mean quality per cycle"))
  lines(after, col = "blue", lwd = 2)
  abline(h = c(20, 30), col = "grey", lty = c(2, 3))
  legend("bottomleft", c("before filtering", "after filtering"),
         col = c("red", "blue"), lwd = 2, bty = "n", cex = 0.8)
}
dev.off()

# base composition per cycle
bases <- c("A", "C", "G", "T", "N")
cols <- c("blue", "orange", "green4", "red", "purple")
open_png("fig2_base_content.png", 1600)
par(mfrow = c(1, 2))
for (r in reads) {
  cc <- rep[[paste0(r, "_before_filtering")]]$content_curves
  plot(100 * cc$A, type = "l", col = cols[1], ylim = c(0, 100),
       xlab = "cycle (bp)", ylab = "base fraction (%)",
       main = paste(r, "base composition per cycle"))
  for (i in 2:5) lines(100 * cc[[bases[i]]], col = cols[i])
  legend("topright", bases, col = cols, lwd = 2, horiz = TRUE, bty = "n", cex = 0.8)
}
dev.off()

# GC content per cycle
open_png("fig3_gc_per_cycle.png")
gc <- 100 * rep$summary$before_filtering$gc_content
plot(100 * rep$read1_before_filtering$content_curves$GC, type = "l", col = "blue",
     ylim = c(0, 100), xlab = "cycle (bp)", ylab = "GC (%)",
     main = paste(sample, "GC content per cycle"))
lines(100 * rep$read2_before_filtering$content_curves$GC, col = "green4", lty = 2)
abline(h = gc, col = "grey", lty = 3)
legend("topright", c("read1", "read2", sprintf("overall %.1f%%", gc)),
       col = c("blue", "green4", "grey"), lty = c(1, 2, 3), bty = "n", cex = 0.8)
dev.off()

# read fate during filtering and base quality before vs after
f <- rep$filtering_result
fate <- c(`passed filter` = f$passed_filter_reads, `low quality` = f$low_quality_reads,
          `too many N` = f$too_many_N_reads, `too short` = f$too_short_reads)
q <- rbind(before = 100 * c(rep$summary$before_filtering$q20_rate, rep$summary$before_filtering$q30_rate),
           after = 100 * c(rep$summary$after_filtering$q20_rate, rep$summary$after_filtering$q30_rate))
colnames(q) <- c(">=Q20 bases", ">=Q30 bases")
open_png("fig4_filtering_and_quality.png", 1600)
par(mfrow = c(1, 2), mar = c(5, 6, 4, 2))
bp <- barplot(fate, col = "steelblue", ylab = "reads (n)", las = 1,
              ylim = c(0, max(fate) * 1.2), main = paste(sample, "read fate during filtering"))
text(bp, fate, labels = sprintf("%d\n%.2f%%", fate, 100 * fate / sum(fate)), pos = 3, cex = 0.8)
bp2 <- barplot(q, beside = TRUE, col = c("red", "blue"), ylim = c(0, 115),
               ylab = "bases (%)", main = paste(sample, "base quality before vs after"))
text(bp2, q, labels = sprintf("%.1f", q), pos = 3, cex = 0.8)
legend("bottomleft", rownames(q), fill = c("red", "blue"), bty = "n", cex = 0.8)
dev.off()

# duplication level
dup <- rep$duplication
open_png("fig5_duplication.png")
k <- max(which(dup$histogram > 0))
counts <- dup$histogram[1:k]
bars <- barplot(counts, names.arg = 1:k, col = "seagreen", border = NA,
        ylim = c(0, max(counts) * 1.15),
        xlab = "times a sequence was observed", ylab = "distinct sequences (n)",
        main = sprintf("%s duplication level (overall rate %.2f%%)", sample, 100 * dup$rate))
text(bars, counts, labels = ifelse(counts > 0, counts, ""), pos = 3, cex = 0.8)
dev.off()

# insert size, or a statement when no pair yields one
ins <- rep$insert_size
open_png("fig6_insert_size.png")
if (sum(ins$histogram) == 0) {
  plot.new()
  title(main = paste(sample, "insert size not determinable"))
  text(0.5, 0.6, "No insert size could be determined\nfor any read pair", cex = 1.3)
  text(0.5, 0.3, sprintf("%d read pairs undetermined\n(read 1 and read 2 do not behave as mates)",
                         ins$unknown), col = "red", cex = 0.9)
} else {
  barplot(ins$histogram, names.arg = seq_along(ins$histogram) - 1, border = NA, space = 0,
          col = "purple", xlab = "insert size (bp)", ylab = "read pairs (n)",
          main = sprintf("%s insert size distribution (peak %d bp)", sample, ins$peak))
}
dev.off()

s <- rep$summary
metrics <- c(sample = sample,
             before_total_reads = s$before_filtering$total_reads,
             before_total_bases = s$before_filtering$total_bases,
             before_q20_rate_percent = 100 * s$before_filtering$q20_rate,
             before_q30_rate_percent = 100 * s$before_filtering$q30_rate,
             before_gc_percent = 100 * s$before_filtering$gc_content,
             before_read1_mean_length_bp = s$before_filtering$read1_mean_length,
             before_read2_mean_length_bp = s$before_filtering$read2_mean_length,
             after_total_reads = s$after_filtering$total_reads,
             after_total_bases = s$after_filtering$total_bases,
             after_q20_rate_percent = 100 * s$after_filtering$q20_rate,
             after_q30_rate_percent = 100 * s$after_filtering$q30_rate,
             after_gc_percent = 100 * s$after_filtering$gc_content,
             after_read1_mean_length_bp = s$after_filtering$read1_mean_length,
             after_read2_mean_length_bp = s$after_filtering$read2_mean_length,
             passed_filter_reads = f$passed_filter_reads,
             low_quality_reads = f$low_quality_reads,
             too_many_N_reads = f$too_many_N_reads,
             too_short_reads = f$too_short_reads,
             duplication_rate_percent = 100 * dup$rate,
             insert_size_undetermined_pairs = ins$unknown)
write.table(data.frame(metric = names(metrics), value = metrics),
            file.path("/results/reports", paste0(sample, "_qc_metrics.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

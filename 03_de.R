# deseq2: robo1 ko vs aavs1

library(tximport)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(jsonlite)

setwd("~/robo1_analysis")
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# save pdf + png
save_plot <- function(p, name, w = 6, h = 5) {
  ggsave(paste0("figures/", name, ".pdf"), p, width = w, height = h)
  ggsave(paste0("figures/", name, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
}

cols <- c(AAVS1 = "#2a78d6", KO = "#eb6834")
dir_cols <- c(up = "#e34948", down = "#2a78d6", ns = "#c9c8c3")
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank()))

# samples
samples <- data.frame(
  sample = c("SS1-KO-1", "SS2-KO-2", "SS3-KO-3", "SS4-AAVS1-1", "SS5-AAVS1-2", "SS6-AAVS1-3"),
  short = c("KO-1", "KO-2", "KO-3", "AAVS1-1", "AAVS1-2", "AAVS1-3"),
  condition = factor(rep(c("KO", "AAVS1"), each = 3), levels = c("AAVS1", "KO"))
)
rownames(samples) <- samples$sample

# gene info from the gencode gtf
tx <- read.delim("ref/tx2gene.tsv", header = FALSE, col.names = c("tx", "gene_id", "symbol", "type", "chr"))
genes <- unique(tx[, c("gene_id", "symbol", "type", "chr")])
rownames(genes) <- genes$gene_id

# salmon -> gene level
files <- file.path("salmon", samples$sample, "quant.sf")
names(files) <- samples$sample
txi <- tximport(files, type = "salmon", tx2gene = tx[, 1:2])
tpm <- txi$abundance

# qc numbers from fastp + salmon
qc <- do.call(rbind, lapply(samples$sample, function(s) {
  fp <- fromJSON(paste0("qc/fastp/", s, ".json"))$summary
  sm <- fromJSON(paste0("salmon/", s, "/aux_info/meta_info.json"))
  lf <- fromJSON(paste0("salmon/", s, "/lib_format_counts.json"))
  data.frame(sample = s,
             read_pairs = fp$before_filtering$total_reads / 2,
             pct_kept = 100 * fp$after_filtering$total_reads / fp$before_filtering$total_reads,
             q30_after = 100 * fp$after_filtering$q30_rate,
             gc = 100 * fp$after_filtering$gc_content,
             pairs_after_trim = sm$num_processed,
             mapped = sm$num_mapped,
             pct_mapped = sm$percent_mapped,
             lib_type = sm$library_types,
             strand_bias = lf$strand_mapping_bias,
             frag_len = sm$frag_length_mean)
}))
qc$condition <- samples$condition
qc$short <- samples$short

# deseq2
dds <- DESeqDataSetFromTximport(txi, samples, ~ condition)

# keep genes with >= 10 reads in at least 3 samples (one group's worth)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)
qc$counts_kept_genes <- colSums(counts(dds))
write.csv(qc, "results/qc_summary.csv", row.names = FALSE)

res <- results(dds, name = "condition_KO_vs_AAVS1", alpha = 0.05)
shr <- lfcShrink(dds, coef = "condition_KO_vs_AAVS1", type = "apeglm")
summary(res)

# full table (log2FC = apeglm shrunk, log2FC_mle = unshrunk)
norm <- counts(dds, normalized = TRUE)
colnames(norm) <- paste0("norm_", samples$short)
gt <- tpm[rownames(dds), ]
colnames(gt) <- paste0("tpm_", samples$short)
de <- data.frame(gene_id = rownames(res), genes[rownames(res), c("symbol", "type", "chr")],
                 baseMean = res$baseMean, log2FC = shr$log2FoldChange, lfcSE = shr$lfcSE,
                 log2FC_mle = res$log2FoldChange, stat = res$stat, pvalue = res$pvalue, padj = res$padj,
                 tpm_AAVS1 = rowMeans(gt[, 4:6]), tpm_KO = rowMeans(gt[, 1:3]),
                 round(norm, 1), round(gt, 2), check.names = FALSE)
de <- de[order(de$padj, de$pvalue), ]
de$direction <- "ns"
de$direction[which(de$padj < 0.05 & de$log2FC > 1)] <- "up"
de$direction[which(de$padj < 0.05 & de$log2FC < -1)] <- "down"

sig <- de[which(de$padj < 0.05), ]
up <- de[de$direction == "up", ]
down <- de[de$direction == "down", ]
write.csv(de, "results/de_all.csv", row.names = FALSE)
write.csv(sig, "results/de_sig_padj05_anyFC.csv", row.names = FALSE)
write.csv(up, "results/de_sig_up.csv", row.names = FALSE)
write.csv(down, "results/de_sig_down.csv", row.names = FALSE)

cat("genes tested:", nrow(de), "\n")
cat("padj < 0.05:", nrow(sig), "(up", sum(sig$log2FC > 0), "/ down", sum(sig$log2FC < 0), ")\n")
cat("padj < 0.05 & |lfc| > 1: up", nrow(up), "/ down", nrow(down), "\n")
cat("protein coding: up", sum(up$type == "protein_coding"), "/ down", sum(down$type == "protein_coding"), "\n")

# qc plots
# library sizes
lib <- rbind(data.frame(short = qc$short, condition = qc$condition, what = "read pairs (raw)", n = qc$read_pairs),
             data.frame(short = qc$short, condition = qc$condition, what = "mapped (salmon)", n = qc$mapped))
lib$what <- factor(lib$what, levels = c("read pairs (raw)", "mapped (salmon)"))
lib$short <- factor(lib$short, levels = samples$short)
p <- ggplot(lib, aes(short, n / 1e6, fill = condition, alpha = what)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = cols) + scale_alpha_manual(values = c(0.45, 1), name = NULL) +
  labs(x = NULL, y = "million fragments", title = "library size") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, "qc_library_sizes", 6, 4)

# mapping rate
qc$short <- factor(qc$short, levels = samples$short)
p <- ggplot(qc, aes(short, pct_mapped, fill = condition)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct_mapped)), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = cols) + coord_cartesian(ylim = c(0, 100)) +
  labs(x = NULL, y = "% fragments mapped", title = "salmon mapping rate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p, "qc_mapping_rate", 5, 4)

# pca on vst
vsd <- vst(dds, blind = TRUE)
pc <- plotPCA(vsd, intgroup = "condition", returnData = TRUE, ntop = 500)
pv <- round(100 * attr(pc, "percentVar"), 1)
pc$short <- samples[rownames(pc), "short"]
p <- ggplot(pc, aes(PC1, PC2, colour = condition)) +
  geom_point(size = 3.5) +
  geom_text_repel(aes(label = short), size = 3.2, show.legend = FALSE) +
  scale_colour_manual(values = cols) +
  labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"), title = "PCA, top 500 variable genes (vst)")
save_plot(p, "qc_pca", 5.5, 4.5)

# sample distance heatmap
d <- as.matrix(dist(t(assay(vsd))))
rownames(d) <- colnames(d) <- samples$short
ann <- data.frame(condition = samples$condition, row.names = samples$short)
hm <- pheatmap(d, annotation_row = ann, annotation_col = ann, annotation_colors = list(condition = cols),
               color = colorRampPalette(c("#104281", "#6da7ec", "#f7f7f5"))(100),
               display_numbers = TRUE, number_format = "%.0f", main = "sample distances (vst, euclidean)", silent = TRUE)
save_plot(hm$gtable, "qc_sample_distances", 6, 5)

# dispersion
disp <- as.data.frame(mcols(dds)[, c("baseMean", "dispGeneEst", "dispFit", "dispersion")])
disp <- disp[disp$baseMean > 0, ]
# gene estimates stuck at the 1e-8 floor squash the plot, leave them out and say how many
low <- sum(disp$dispGeneEst < 1e-6)
disp$dispGeneEst[disp$dispGeneEst < 1e-6] <- NA
p <- ggplot(disp, aes(baseMean)) +
  geom_point(aes(y = dispGeneEst, colour = "gene estimate"), size = 0.3, alpha = 0.4, na.rm = TRUE) +
  geom_point(aes(y = dispersion, colour = "final (shrunk)"), size = 0.3, alpha = 0.4) +
  geom_line(aes(y = dispFit, colour = "fitted trend"), linewidth = 0.8) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c("gene estimate" = "grey55", "final (shrunk)" = "#2a78d6", "fitted trend" = "#e34948"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "mean normalized count", y = "dispersion", title = "DESeq2 dispersion estimates",
       caption = paste(low, "genes with gene-wise estimate at the lower bound not shown"))
save_plot(p, "qc_dispersion", 6, 4.5)

# ma plot
p <- ggplot(de, aes(baseMean, log2FC, colour = direction)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
  scale_x_log10() + scale_colour_manual(values = dir_cols) +
  labs(x = "mean normalized count", y = "log2 fold change (KO / AAVS1, shrunk)", title = "MA plot")
save_plot(p, "ma_plot", 6, 4.5)

# volcano, label top 15 each way + robo1
de$nlp <- -log10(de$padj)
lab <- c(head(up$gene_id, 15), head(down$gene_id, 15), de$gene_id[de$symbol == "ROBO1"])
de$label <- ifelse(de$gene_id %in% lab, de$symbol, NA)
p <- ggplot(de[!is.na(de$padj), ], aes(log2FC, nlp, colour = direction)) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  geom_text_repel(aes(label = label), colour = "grey10", size = 2.8, max.overlaps = 50,
                  min.segment.length = 0, segment.colour = "grey50", na.rm = TRUE) +
  scale_colour_manual(values = dir_cols,
                      labels = c(down = paste0("down (", nrow(down), ")"), ns = "ns", up = paste0("up (", nrow(up), ")"))) +
  scale_y_sqrt(breaks = c(1.3, 5, 20, 50, 100, 200), labels = c("1.3", "5", "20", "50", "100", "200")) +
  labs(x = "log2 fold change (KO / AAVS1, shrunk)", y = "-log10 padj (sqrt scale)", colour = NULL,
       title = "ROBO1 KO vs AAVS1", subtitle = "padj < 0.05 and |log2FC| > 1")
save_plot(p, "volcano", 7, 6)

# heatmap of top de genes (up to 25 up + 25 down by padj), row z-score of vst
top <- c(head(up$gene_id, 25), head(down$gene_id, 25))
ttl <- paste0("top DE genes: ", min(25, nrow(up)), " up + ", min(25, nrow(down)), " down (row z-score, vst)")
m <- assay(vsd)[top, ]
m <- t(scale(t(m)))
rownames(m) <- genes[top, "symbol"]
colnames(m) <- samples$short
hm <- pheatmap(m, annotation_col = ann, annotation_colors = list(condition = cols),
               color = colorRampPalette(c("#1c5cab", "#f0efec", "#c23b3a"))(100), breaks = seq(-2, 2, length.out = 101),
               cluster_cols = TRUE, fontsize_row = 7, main = ttl, silent = TRUE)
save_plot(hm$gtable, "heatmap_top_de", 5.5, 9)

# marker checks
mk <- data.frame(
  symbol = c("ROBO1", "ROBO2", "ROBO3", "ROBO4", "SLIT1", "SLIT2", "SLIT3",
             "PPP1R12C",
             "SOX2", "OLIG2", "NES", "PROM1", "BMI1", "NANOG", "POU3F2", "SALL2",
             "MKI67", "TOP2A",
             "XIST", "RPS4Y1", "DDX3Y"),
  group = c(rep("slit/robo", 7), "aavs1 locus", rep("stemness", 8), rep("proliferation", 2), rep("sex check", 3))
)
ids <- genes$gene_id[match(mk$symbol, genes$symbol)]
mk <- cbind(mk, round(tpm[ids, ], 2))
colnames(mk)[3:8] <- paste0("tpm_", samples$short)
mk$tpm_AAVS1 <- rowMeans(mk[, 6:8])
mk$tpm_KO <- rowMeans(mk[, 3:5])
mk$log2FC <- de$log2FC[match(ids, de$gene_id)]
mk$padj <- de$padj[match(ids, de$gene_id)]
mk$note <- ifelse(ids %in% de$gene_id, "", "filtered (low counts)")
write.csv(mk, "results/marker_check.csv", row.names = FALSE)

# long format for plotting
mk_long <- function(g) {
  x <- mk[mk$symbol %in% g, ]
  out <- data.frame(symbol = rep(x$symbol, 6), short = rep(samples$short, each = nrow(x)),
                    tpm = unlist(x[, 3:8]), condition = rep(samples$condition, each = nrow(x)))
  stat <- ifelse(is.na(x$padj), "not tested", sprintf("LFC %.2f, padj %.2g", x$log2FC, x$padj))
  out$panel <- factor(paste0(out$symbol, "\n", stat[match(out$symbol, x$symbol)]),
                      levels = paste0(x$symbol, "\n", stat)[match(g, x$symbol)])
  out
}

mplot <- function(d, title) {
  ggplot(d, aes(condition, tpm, colour = condition)) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3, colour = "grey30") +
    geom_point(size = 2.2, position = position_jitter(width = 0.08, seed = 1)) +
    facet_wrap(~ panel, scales = "free_y", nrow = 2) +
    scale_colour_manual(values = cols) + expand_limits(y = 0) +
    labs(x = NULL, y = "TPM", title = title) +
    theme(legend.position = "none", strip.text = element_text(size = 8))
}

p <- mplot(mk_long(c("SOX2", "OLIG2", "NES", "PROM1", "BMI1", "NANOG", "POU3F2", "SALL2", "MKI67", "TOP2A")),
           "stemness + proliferation markers")
save_plot(p, "stemness_markers", 10, 5)

p <- mplot(mk_long(c("ROBO1", "ROBO2", "ROBO3", "ROBO4", "SLIT1", "SLIT2", "SLIT3", "PPP1R12C")),
           "SLIT/ROBO family + AAVS1 locus (PPP1R12C)")
save_plot(p, "robo_slit_markers", 9, 5)

# robo1 transcripts, to see if any isoform drops (nmd) even if gene total doesn't
txo <- tximport(files, type = "salmon", txOut = TRUE)
rtx <- tx[tx$symbol == "ROBO1", "tx"]
rt <- round(txo$abundance[rtx, ], 2)
colnames(rt) <- paste0("tpm_", samples$short)
rt <- data.frame(tx = rtx, rt, check.names = FALSE)
rt <- rt[rowSums(rt[, -1]) > 0, ]
rt$tpm_AAVS1 <- rowMeans(rt[, 5:7])
rt$tpm_KO <- rowMeans(rt[, 2:4])
rt <- rt[order(-rt$tpm_AAVS1), ]
write.csv(rt, "results/robo1_transcripts.csv", row.names = FALSE)

writeLines(capture.output(sessionInfo()), "results/sessionInfo_de.txt")

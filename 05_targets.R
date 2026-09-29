# car-t co-target candidates (surface genes up in ko) + secreted de genes + excel workbook

library(ggplot2)
library(ggrepel)
library(httr)
library(openxlsx)

setwd("~/robo1_analysis")
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank()))

save_plot <- function(p, name, w = 6, h = 5) {
  ggsave(paste0("figures/", name, ".pdf"), p, width = w, height = h)
  ggsave(paste0("figures/", name, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
}

# annotation downloads (hpa v25, gtex v10, surfy/cspa)
dir.create("ref/annot", showWarnings = FALSE)
get <- function(url, f) if (!file.exists(f)) download.file(url, f, mode = "wb", quiet = TRUE)
get("https://www.proteinatlas.org/download/proteinatlas.tsv.zip", "ref/annot/proteinatlas.tsv.zip")
get("https://www.proteinatlas.org/download/tsv/normal_ihc_data.tsv.zip", "ref/annot/normal_ihc_data.tsv.zip")
get("https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz", "ref/annot/gtex_v10_median_tpm.gct.gz")
get("https://wollscheidlab.org/SURFY/table_S3_surfaceome.xlsx", "ref/annot/surfy_table_S3.xlsx")

de <- read.csv("results/de_all.csv", check.names = FALSE)
de$ensg <- sub("\\..*", "", de$gene_id)

# hpa: subcellular location, secretome, protein class
hpa <- read.delim(unz("ref/annot/proteinatlas.tsv.zip", "proteinatlas.tsv"), check.names = FALSE)
hpa <- hpa[match(de$ensg, hpa$Ensembl), ]
de$hpa_location <- hpa$`Subcellular location`
de$hpa_secretome <- hpa$`Secretome location`
de$hpa_secretome_function <- hpa$`Secretome function`
de$hpa_protein_class <- hpa$`Protein class`
de$hpa_tissue_specificity <- hpa$`RNA tissue specificity`
de$description <- hpa$`Gene description`

# surfy (in silico surfaceome) + cspa (mass spec on real cell surfaces), same table
sf <- read.xlsx("ref/annot/surfy_table_S3.xlsx", sheet = 1, startRow = 2)
sf <- sf[!is.na(sf$Ensembl.gene), ]
sf <- sf[match(de$ensg, sf$Ensembl.gene), ]
de$surfy <- sf$Surfaceome.Label
de$cspa <- sf$CSPA.category
de$cd <- sf$CD.number

# surface score: one point each for hpa plasma membrane, surfy surface, cspa high confidence
de$hpa_pm <- grepl("Plasma membrane", de$hpa_location)
de$surface_score <- de$hpa_pm + (de$surfy %in% "surface") + grepl("1 - high", de$cspa)

# gtex normal tissue: max median tpm across brain regions and across everything else (no cell lines)
gt <- read.delim(gzfile("ref/annot/gtex_v10_median_tpm.gct.gz"), skip = 2, check.names = FALSE)
gt <- gt[match(de$ensg, sub("\\..*", "", gt$Name)), ]
tis <- setdiff(colnames(gt), c("Name", "Description"))
brain <- grep("^Brain", tis, value = TRUE)
other <- setdiff(tis, c(brain, grep("^Cells", tis, value = TRUE)))
b <- as.matrix(gt[, brain])
o <- as.matrix(gt[, other])
de$brain_max_tpm <- round(apply(b, 1, max), 2)
de$brain_max_region <- ifelse(is.na(de$brain_max_tpm), NA, brain[max.col(replace(b, is.na(b), 0), "first")])
de$other_max_tpm <- round(apply(o, 1, max), 2)
de$other_max_tissue <- ifelse(is.na(de$other_max_tpm), NA, other[max.col(replace(o, is.na(o), 0), "first")])

# hpa ihc (protein) in normal brain, highest level seen
ihc <- read.delim(unz("ref/annot/normal_ihc_data.tsv.zip", "normal_ihc_data.tsv"), check.names = FALSE)
ihc <- ihc[ihc$Tissue %in% c("Cerebral cortex", "Cerebellum", "Hippocampus", "Caudate", "Hypothalamus",
                             "Dorsal raphe", "Substantia nigra", "Choroid plexus") & ihc$Reliability != "Uncertain", ]
lv <- c("Not detected", "Low", "Medium", "High")
ihc$lvl <- match(ihc$Level, lv)
ihc <- ihc[!is.na(ihc$lvl), ]
bmax <- tapply(ihc$lvl, ihc$Gene, max, na.rm = TRUE)
de$brain_ihc_max <- lv[bmax[de$ensg]]

# candidates: sig up (padj < 0.05) with at least one surface call
# only 12 genes pass the main lfc > 1 cutoff, so smaller lfc are kept too and split into tiers:
# tier 1 = lfc > 1 and >= 2 surface sources, tier 2 = any lfc and >= 2 sources, tier 3 = one source only
# (one source alone is often noise, e.g. hpa IF plasma membrane on kinesins or cspa picking up secreted ecm)
ct <- de[which(de$padj < 0.05 & de$log2FC > 0 & de$surface_score >= 1), ]
ct$tier <- ifelse(ct$surface_score < 2, 3, ifelse(ct$log2FC > 1, 1, 2))

# existing drugs / antibodies from dgidb
# only interactions with a stated type (drops biomarker-only links, e.g. alk <-> nivolumab)
q <- paste0('{ genes(names: [', paste0('"', ct$symbol, '"', collapse = ","), ']) { nodes { name interactions { drug { name approved } interactionTypes { type } } } } }')
r <- POST("https://dgidb.org/api/graphql", body = list(query = q), encode = "json")
nodes <- content(r, as = "parsed")$data$genes$nodes
dg <- do.call(rbind, lapply(nodes, function(n) {
  ints <- Filter(function(i) length(i$interactionTypes) > 0, n$interactions)
  if (length(ints) == 0) return(data.frame(symbol = n$name, n_drugs = 0, approved_drugs = "", antibody_drugs = ""))
  d <- sapply(ints, function(i) i$drug$name)
  ok <- sapply(ints, function(i) isTRUE(i$drug$approved))
  ab <- grepl("MAB\\b", d) | sapply(ints, function(i) "antibody" %in% sapply(i$interactionTypes, function(t) t$type))
  data.frame(symbol = n$name, n_drugs = length(unique(d)),
             approved_drugs = paste(head(unique(d[ok]), 5), collapse = "; "),
             antibody_drugs = paste(head(unique(d[ab]), 5), collapse = "; "))
}))
m <- match(ct$symbol, dg$symbol)
ct$n_drugs <- ifelse(is.na(m), 0, dg$n_drugs[m])
ct$approved_drugs <- ifelse(is.na(m), "", dg$approved_drugs[m])
ct$antibody_drugs <- ifelse(is.na(m), "", dg$antibody_drugs[m])

# surface targets with car-t / bispecific / adc programs (hand list from the literature)
known <- c("EGFR", "IL13RA2", "ERBB2", "EPHA2", "CD276", "CD70", "PROM1", "MMP2", "CSPG4", "PDPN", "GPC2", "GPC3",
           "ROBO1", "MSLN", "CLDN18", "ROR1", "MUC1", "CEACAM5", "L1CAM", "DLL3", "FOLH1", "PSCA", "FAP", "EPCAM",
           "TACSTD2", "NECTIN4", "FOLR1", "MICA", "MICB", "ULBP2", "CD44", "AXL", "MET")
ct$known_immunotherapy_target <- ct$symbol %in% known
ct$fda_drug_target <- grepl("FDA approved drug targets", ct$hpa_protein_class)

# rank: average percentile of big lfc, strong padj, high ko expression, surface score, low brain, low other tissue
pr <- function(x) rank(x, na.last = "keep") / sum(!is.na(x))
ct$score <- round(rowMeans(cbind(pr(ct$log2FC), pr(-log10(ct$padj)), pr(ct$tpm_KO), ct$surface_score / 3,
                                 pr(-ct$brain_max_tpm), pr(-ct$other_max_tpm)), na.rm = TRUE), 3)
ct <- ct[order(ct$tier, -ct$score), ]
ct$rank <- seq_len(nrow(ct))

cols <- c("rank", "tier", "symbol", "description", "score", "log2FC", "padj", "tpm_AAVS1", "tpm_KO",
          "surface_score", "hpa_pm", "surfy", "cspa", "cd", "hpa_location",
          "brain_max_tpm", "brain_max_region", "brain_ihc_max", "other_max_tpm", "other_max_tissue", "hpa_tissue_specificity",
          "hpa_secretome", "known_immunotherapy_target", "fda_drug_target", "n_drugs", "approved_drugs", "antibody_drugs", "gene_id", "type")
ct <- ct[, cols]
write.csv(ct, "results/cotarget_candidates.csv", row.names = FALSE)
cat("surface candidates:", nrow(ct), "| tier 1", sum(ct$tier == 1), "/ tier 2", sum(ct$tier == 2), "/ tier 3", sum(ct$tier == 3),
    "| main cutoff (lfc > 1):", sum(ct$log2FC > 1), "\n")
print(head(ct[, c("rank", "tier", "symbol", "log2FC", "padj", "tpm_KO", "surface_score", "brain_max_tpm", "other_max_tpm")], 15))

# co-target plot: lfc vs normal brain expression
ct$lab <- ifelse(ct$rank <= 15, ct$symbol, NA)
ct$surface <- factor(ct$surface_score, levels = 1:3)
p <- ggplot(ct, aes(log10(brain_max_tpm + 1), log2FC)) +
  geom_vline(xintercept = log10(11), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(aes(size = tpm_KO, colour = surface), alpha = 0.85) +
  geom_text_repel(aes(label = lab), size = 3, max.overlaps = 30, min.segment.length = 0, na.rm = TRUE) +
  scale_colour_manual(values = c("1" = "#86b6ef", "2" = "#2a78d6", "3" = "#104281"), name = "surface\nevidence (0-3)", drop = FALSE) +
  scale_size_continuous(range = c(1.2, 7), name = "TPM in KO", trans = "sqrt") +
  scale_x_continuous(breaks = log10(c(0, 1, 10, 100, 1000) + 1), labels = c(0, 1, 10, 100, 1000)) +
  labs(x = "max median TPM in normal brain (GTEx)", y = "log2 fold change (KO / AAVS1)",
       title = "co-target candidates: surface genes up in ROBO1 KO", subtitle = "all padj < 0.05; top 15 (by tier, then score) labelled; lines at log2FC = 1 and 10 TPM")
save_plot(p, "cotarget_scatter", 8, 6)

# top 20 profile, one panel per measure
t20 <- head(ct, 20)
lv20 <- rev(t20$symbol)
long <- rbind(
  data.frame(symbol = t20$symbol, what = "log2FC (KO/AAVS1)", value = t20$log2FC),
  data.frame(symbol = t20$symbol, what = "TPM in KO", value = t20$tpm_KO),
  data.frame(symbol = t20$symbol, what = "max TPM normal brain", value = t20$brain_max_tpm),
  data.frame(symbol = t20$symbol, what = "max TPM other tissue", value = t20$other_max_tpm)
)
long$symbol <- factor(long$symbol, levels = lv20)
long$what <- factor(long$what, levels = unique(long$what))
p <- ggplot(long, aes(value, symbol)) +
  geom_segment(aes(x = 0, xend = value, yend = symbol), colour = "grey70") +
  geom_point(size = 2.5, colour = "#2a78d6") +
  facet_wrap(~ what, nrow = 1, scales = "free_x") +
  labs(x = NULL, y = NULL, title = "top 20 co-target candidates", subtitle = "ordered by tier, then combined score")
save_plot(p, "cotarget_top20", 11, 6)

# secreted de genes (hpa secretome), cells were grown alone so these are what the tumour would put out
sec <- de[de$direction != "ns" & grepl("^Secreted", de$hpa_secretome), ]
sec <- sec[order(-sec$log2FC), c("symbol", "description", "direction", "log2FC", "padj", "tpm_AAVS1", "tpm_KO",
                                 "hpa_secretome", "hpa_secretome_function", "hpa_protein_class", "brain_max_tpm", "gene_id", "type")]
write.csv(sec, "results/secreted_de.csv", row.names = FALSE)
cat("secreted de genes: up", sum(sec$direction == "up"), "/ down", sum(sec$direction == "down"), "\n")

s <- rbind(head(sec[sec$direction == "up", ], 15), tail(sec[sec$direction == "down", ], 15))
s$symbol <- factor(s$symbol, levels = rev(unique(s$symbol)))
p <- ggplot(s, aes(log2FC, symbol, colour = direction)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_segment(aes(x = 0, xend = log2FC, yend = symbol)) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c(up = "#e34948", down = "#2a78d6")) +
  labs(x = "log2 fold change (KO / AAVS1)", y = NULL, colour = NULL, title = "secreted DE genes (HPA secretome)",
       subtitle = paste0("padj < 0.05, |log2FC| > 1: ", sum(s$direction == "up"), " up, ", sum(s$direction == "down"), " down (max 15 each shown)"))
save_plot(p, "secreted_de", 6.5, 7)

# established car-t / antibody targets: are they still on the ko cells? (a co-target only has to stay present)
kt <- de[de$symbol %in% known, c("symbol", "log2FC", "padj", "tpm_AAVS1", "tpm_KO", "surface_score",
                                 "brain_max_tpm", "other_max_tpm", "other_max_tissue", "gene_id")]
kt <- kt[order(-kt$tpm_KO), ]
write.csv(kt, "results/known_targets_check.csv", row.names = FALSE)

kt$lab <- paste0(kt$symbol, ifelse(!is.na(kt$padj) & kt$padj < 0.05, " *", ""))
kl <- rbind(data.frame(lab = kt$lab, condition = "AAVS1", tpm = kt$tpm_AAVS1),
            data.frame(lab = kt$lab, condition = "KO", tpm = kt$tpm_KO))
kl$lab <- factor(kl$lab, levels = rev(kt$lab))
p <- ggplot(kl, aes(tpm + 1, lab)) +
  geom_line(aes(group = lab), colour = "grey70") +
  geom_point(aes(colour = condition), size = 2.5) +
  scale_x_log10() +
  scale_colour_manual(values = c(AAVS1 = "#2a78d6", KO = "#eb6834"), name = NULL) +
  labs(x = "TPM + 1 (log scale)", y = NULL, title = "established immunotherapy targets: AAVS1 vs ROBO1 KO",
       subtitle = "* padj < 0.05; genes not expressed in these cells are left out")
save_plot(p, "known_targets", 6.5, 6)

# one excel workbook, a sheet per table
wb <- createWorkbook()
tabs <- c(qc_summary = "qc_summary", de_all = "de_all", sig_up = "de_sig_up", sig_down = "de_sig_down",
          sig_padj05_anyFC = "de_sig_padj05_anyFC", markers = "marker_check", robo1_transcripts = "robo1_transcripts",
          gsea_hallmark = "gsea_hallmark", gsea_go_bp = "gsea_go_bp", gsea_reactome = "gsea_reactome",
          gsea_kegg = "gsea_kegg", gsea_gbm = "gsea_gbm", cotargets = "cotarget_candidates", known_targets = "known_targets_check",
          secreted = "secreted_de")
for (n in names(tabs)) {
  addWorksheet(wb, n)
  writeData(wb, n, read.csv(paste0("results/", tabs[n], ".csv"), check.names = FALSE))
  freezePane(wb, n, firstRow = TRUE)
}
saveWorkbook(wb, "results/robo1_ko_rnaseq.xlsx", overwrite = TRUE)

writeLines(capture.output(sessionInfo()), "results/sessionInfo_targets.txt")

# gsea with fgsea: shrunk lfc ranked list (planned) + wald stat ranked list

library(fgsea)
library(msigdbr)
library(openxlsx)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)

setwd("~/robo1_analysis")
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank()))

save_plot <- function(p, name, w = 6, h = 5) {
  ggsave(paste0("figures/", name, ".pdf"), p, width = w, height = h)
  ggsave(paste0("figures/", name, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
}

de <- read.csv("results/de_all.csv", check.names = FALSE)

# ranked lists, one value per symbol (keep the copy with the biggest |value|)
rank_by <- function(col) {
  d <- de[!is.na(de[[col]]) & de$symbol != "", ]
  d <- d[order(-abs(d[[col]])), ]
  d <- d[!duplicated(d$symbol), ]
  sort(setNames(d[[col]], d$symbol), decreasing = TRUE)
}
# apeglm pulls ~3/4 of genes to lfc ~0, so the shrunk lfc list is mostly ties and gsea has little power.
# the wald stat (usual fgsea input) keeps the ordering of all genes, so run both and report both
ranks <- list(wald_stat = rank_by("stat"), shrunk_lfc = rank_by("log2FC"))
cat("genes with |shrunk lfc| < 0.05:", sum(abs(ranks$shrunk_lfc) < 0.05), "of", length(ranks$shrunk_lfc), "\n")

# msigdb sets
msig <- function(coll, sub = NULL) {
  m <- msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub)
  split(m$gene_symbol, m$gs_name)
}
sets <- list(
  hallmark = msig("H"),
  go_bp = msig("C5", "GO:BP"),
  reactome = msig("C2", "CP:REACTOME"),
  kegg = msig("C2", "CP:KEGG_LEGACY")
)

# gbm sets: neftel 2019 cell states (table s2), verhaak 2010 subtypes, 3ca glioma metaprograms, stemness sets
dir.create("ref/annot", showWarnings = FALSE)
nf <- "ref/annot/neftel2019_tableS2.xlsx"
if (!file.exists(nf)) download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867419306877-mmc2.xlsx", nf, mode = "wb", quiet = TRUE)
neftel <- read.xlsx(nf, startRow = 4)
neftel <- lapply(neftel, function(x) x[!is.na(x)])
names(neftel) <- paste0("NEFTEL_", gsub("[/.]", "", names(neftel)))

cgp <- c(msig("C2", "CGP"), msig("C4", "3CA"))
keep <- c("VERHAAK_GLIOBLASTOMA_CLASSICAL", "VERHAAK_GLIOBLASTOMA_MESENCHYMAL", "VERHAAK_GLIOBLASTOMA_NEURAL", "VERHAAK_GLIOBLASTOMA_PRONEURAL",
          "GAVISH_3CA_MALIGNANT_METAPROGRAM_16_MES_GLIOMA", "GAVISH_3CA_MALIGNANT_METAPROGRAM_26_NPC_GLIOMA",
          "GAVISH_3CA_MALIGNANT_METAPROGRAM_25_ASTROCYTES", "GAVISH_3CA_MALIGNANT_METAPROGRAM_27_OLIGO_PROGENITOR",
          "GAVISH_3CA_MALIGNANT_METAPROGRAM_29_NPC_OPC",
          "BEIER_GLIOMA_STEM_CELL_UP", "BEIER_GLIOMA_STEM_CELL_DN", "ZHENG_GLIOBLASTOMA_PLASTICITY_UP", "ZHENG_GLIOBLASTOMA_PLASTICITY_DN",
          "BENPORATH_ES_1", "BENPORATH_ES_2", "BENPORATH_SOX2_TARGETS", "BENPORATH_NANOG_TARGETS", "BENPORATH_OCT4_TARGETS",
          "WONG_EMBRYONIC_STEM_CELL_CORE", "MALTA_CURATED_STEMNESS_MARKERS", "RAMALHO_STEMNESS_UP", "RAMALHO_STEMNESS_DN")
gbm <- c(neftel, cgp[keep])

# older symbols (e.g. neftel 2019) -> current symbols
fix_sym <- function(g) {
  miss <- !(g %in% names(ranks$wald_stat))
  new <- suppressMessages(mapIds(org.Hs.eg.db, g[miss], "SYMBOL", "ALIAS", multiVals = "first"))
  g[miss] <- ifelse(is.na(new), g[miss], new)
  unique(g)
}
sets$gbm <- lapply(gbm, fix_sym)

# run fgsea
run <- function(s, r, min = 15) {
  set.seed(1)
  x <- fgsea(s, r, minSize = min, maxSize = 500, eps = 0, nproc = 8)
  x[order(x$pval), ]
}
# main = not explained by a more significant overlapping term (for go / reactome)
collapse <- function(x, s, r) {
  sig <- x[x$padj < 0.05, ]
  if (nrow(sig) == 0) return(rep(FALSE, nrow(x)))
  x$pathway %in% collapsePathways(sig, s, r)$mainPathways
}

gs <- list()
for (n in names(sets)) {
  res <- list()
  for (k in names(ranks)) {
    x <- run(sets[[n]], ranks[[k]], ifelse(n == "gbm", 10, 15))
    x$ranking <- k
    x$main <- if (n %in% c("go_bp", "reactome")) collapse(x, sets[[n]], ranks[[k]]) else NA
    res[[k]] <- x
  }
  gs[[n]] <- rbind(res$wald_stat, res$shrunk_lfc)
}

# save tables
for (n in names(gs)) {
  r <- as.data.frame(gs[[n]])
  r$leadingEdge <- sapply(r$leadingEdge, paste, collapse = ";")
  write.csv(r, paste0("results/gsea_", n, ".csv"), row.names = FALSE)
  for (k in names(ranks)) {
    x <- r[r$ranking == k, ]
    cat(n, k, ": padj < 0.05 up", sum(x$padj < 0.05 & x$NES > 0), "/ down", sum(x$padj < 0.05 & x$NES < 0), "\n")
  }
}

# dot plots: top 10 each direction by the wald ranking (plus anything sig by shrunk lfc), both rankings side by side
nice <- function(x) {
  x <- sub("^GAVISH_3CA_MALIGNANT_METAPROGRAM_[0-9]+_", "3CA_", x)
  x <- sub("^(HALLMARK|GOBP|REACTOME|KEGG)_", "", x)
  x <- tolower(gsub("_", " ", x))
  ifelse(nchar(x) > 55, paste0(substr(x, 1, 52), "..."), x)
}
rank_lab <- c(wald_stat = "ranked by Wald stat", shrunk_lfc = "ranked by shrunk LFC (planned)")
dir_cols <- c("up in KO" = "#e34948", "down in KO" = "#2a78d6")

dotplot <- function(r, title, n = 10) {
  r <- as.data.frame(r)
  w <- r[r$ranking == "wald_stat" & r$padj < 0.05 & (is.na(r$main) | r$main), ]
  l <- r[r$ranking == "shrunk_lfc" & r$padj < 0.05, ]
  pick <- unique(c(head(w$pathway[w$NES > 0], n), head(w$pathway[w$NES < 0], n), head(l$pathway, 5)))
  d <- r[r$pathway %in% pick, ]
  d$name <- nice(d$pathway)
  o <- d[d$ranking == "wald_stat", ]
  d$name <- factor(d$name, levels = unique(o$name[order(o$NES)]))
  d$ranking <- factor(rank_lab[d$ranking], levels = rank_lab)
  d$dir <- ifelse(d$NES > 0, "up in KO", "down in KO")
  d$sig <- ifelse(d$padj < 0.05, "padj < 0.05", "ns")
  ggplot(d, aes(NES, name, size = -log10(padj), colour = dir, shape = sig)) +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_point(stroke = 0.8) +
    facet_wrap(~ ranking) +
    scale_colour_manual(values = dir_cols, name = NULL) +
    scale_shape_manual(values = c("padj < 0.05" = 16, "ns" = 1), name = NULL) +
    scale_size_continuous(range = c(1.5, 6), name = "-log10 padj") +
    labs(x = "normalized enrichment score", y = NULL, title = title)
}
titles <- c(hallmark = "Hallmark", go_bp = "GO biological process (collapsed)", reactome = "Reactome (collapsed)", kegg = "KEGG (legacy)")
for (n in names(titles)) {
  p <- dotplot(gs[[n]], paste0("GSEA: ", titles[n]))
  save_plot(p, paste0("gsea_", n), 11, 7)
}

# gbm sets: all of them, significant or not
g <- as.data.frame(gs$gbm)
g$name <- nice(g$pathway)
o <- g[g$ranking == "wald_stat", ]
g$name <- factor(g$name, levels = o$name[order(o$NES)])
g$ranking <- factor(rank_lab[g$ranking], levels = rank_lab)
g$dir <- ifelse(g$NES > 0, "up in KO", "down in KO")
g$sig <- ifelse(g$padj < 0.05, "padj < 0.05", "ns")
p <- ggplot(g, aes(NES, name, fill = dir, alpha = sig)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, colour = "grey40") +
  facet_wrap(~ ranking) +
  scale_fill_manual(values = dir_cols, name = NULL) +
  scale_alpha_manual(values = c("padj < 0.05" = 1, "ns" = 0.3), name = NULL) +
  labs(x = "normalized enrichment score", y = NULL, title = "GBM cell states, subtypes and stemness signatures")
save_plot(p, "gsea_gbm", 11, 7)

# running-sum plots for the top gbm sets (wald ranking)
w <- gs$gbm[gs$gbm$ranking == "wald_stat", ]
top <- head(w$pathway, 8)
pw <- setNames(sets$gbm[top], nice(top))
w$pathway <- nice(w$pathway)
p <- plotGseaTable(pw, ranks$wald_stat, w, gseaParam = 0.5, colwidths = c(5, 3, 0.8, 1.2, 1.2))
save_plot(p, "gsea_gbm_table", 10, 5)

writeLines(capture.output(sessionInfo()), "results/sessionInfo_pathways.txt")

# ROBO1 KO vs AAVS1 (BT241) results notes

## QC
- 15.3-18.9M read pairs/sample (51bp PE). fastp kept 99.3% of reads, Q30 ~94.6%, adapters <0.3%, dup 15-18%, GC 50.6%. nothing weird between samples
- salmon mapping 95.2-95.7%, lib type ISR in all 6 (dUTP stranded), frag length ~220bp
- PCA: PC1 51.5%, KO and AAVS1 separate completely on PC1. PC2 (15%) is just spread within groups. sample distance heatmap also clusters by group (within 10-11, between 12-13)
- XIST expressed, no Y genes, so female and the same in all 6
- PPP1R12C (AAVS1 locus) 58 vs 60 TPM, not affected by the control cut

## is the KO real
- ROBO1 mRNA down: log2FC -2.49, padj 6e-35, 15.8 -> 2.8 TPM (~82% lower). all expressed isoforms go down, probably NMD. so the KO shows at RNA level as well as on the western
- ROBO2 also down (log2FC -2.27, padj 7e-18, 4.3 -> 0.9 TPM). ROBO2 is ~1Mb from ROBO1 on chr3p12, need to figure out why (see caveats)
- ROBO3, ROBO4 no change. SLIT1-3 barely expressed (<1 TPM)

## DE
16,912 genes tested after filtering
- padj < 0.05: 803 (409 up / 394 down)
- padj < 0.05 and |log2FC| > 1: 12 up / 68 down (8 / 65 protein coding)
- mostly small changes, and the big changes are nearly all decreases
- top down by padj: CHI3L1, IGFBP5, SERPINE2, BST2, NEFL, PCDH19, ROBO1, THY1, ROBO2, NNMT
- biggest fold decreases (>20x): SLITRK2, CDH11, CPNE4, SPON1, BST2, POSTN
- the 12 up: ANPEP, RENO1, PLCXD3, RIMS1, ADAMTS6, EGF, KIF17, MYH9-DT + 4 unnamed/readthrough loci at low TPM (ENSG00000310789, ENSG00000268790, ABCF2-H2BK1, ENSG00000275180). not trusting those 4

## stemness + proliferation markers
- no change in SOX2, OLIG2, NES, BMI1, SALL2. PROM1 and NANOG are barely expressed anyway
- POU3F2 down (log2FC -0.86, padj 4e-4) but low expression (2.8 -> 1.5 TPM)
- MKI67, TOP2A flat. proliferation only shows up at pathway level: Neftel G2M down (padj 0.004), hallmark G2M slightly down (padj 0.04)

## GSEA
- ranking by apeglm shrunk LFC (original plan) gave almost nothing: hallmark IFN-a, IFN-g and OXPHOS down, plus 1 GBM set. apeglm shrinks ~73% of genes to ~0 so most of the list is ties
- reran with the DESeq2 Wald stat as the ranking. sig sets (padj < 0.05, up/down): hallmark 13/7, GO BP 118/72, reactome 57/23, KEGG 13/6, GBM sets 5/11. whatever the LFC ranking did find matches the Wald results
- down in KO:
  - type I/II IFN response (IFN-a NES -2.69, padj 6e-12). leading edge BST2, IFI27, IFITM1/3, OAS1, RSAD2, ISG15, HLA-A/B/C a bit
  - OXPHOS / resp chain (NES -2.43, padj 6e-12)
  - ECM / EMT (SERPINE2, THY1, TNC, FN1, CDH11)
  - GBM programs: Neftel NPC1, OPC, AC, MES1, G2M. 3CA OPC, NPC, astro. Wong ES core. Verhaak classical
- up in KO:
  - TNFa/NF-kB (NES 2.14, padj 6e-8), TGF-b (2.13), hypoxia (1.97), apoptosis, KRAS
  - RTK / MAPK / ErbB / PI3K-AKT and focal adhesion (reactome, KEGG)
  - Neftel MES2, Verhaak mesenchymal, 3CA MES glioma, Beier glioma stem cell DN (genes that are normally low in GSCs)
- my interpretation: without ROBO1 the cells lose the developmental/progenitor-like programs (NPC/OPC/AC) and move toward a stress/injury-like MES2 state, with less IFN and OXPHOS. fits with fewer secondary spheres even though SOX2/OLIG2/NES don't change
- SLIT/ROBO reactome sets are not enriched

## co-targets (results/cotarget_candidates.csv)
- 139 genes up (padj < 0.05) with at least 1 surface annotation
- only 3 of these have log2FC > 1, and ANPEP is the only one with good surface evidence
- tiers: 1 = ANPEP, 2 = 47 genes (2+ surface sources, smaller FC), 3 = 91 genes (1 source only, mostly not real surface proteins)
- top 10:
  1. ANPEP (CD13): 2.9 -> 8.3 TPM, brain 1.5 TPM, but ~1000 TPM in gut/kidney
  2. NTNG1: 9.9 -> 19.5, brain 5.5, low in other tissues
  3. NRP2: 163 -> 229, IHC high in brain
  4. NECTIN3 (CD113)
  5. SLC14A1
  6. ALK (CD246), ALK inhibitors available
  7. CDH2 (N-cadherin), high in brain
  8. TNFRSF10D (CD264)
  9. ANO1
  10. PCDH1
- #11 IL13RA2 (already a GBM CAR-T target, 38 -> 50 TPM, padj 0.002), #12 CD274/PD-L1 (10 -> 14 TPM)
- most of these fold changes are small (1.2-1.9x)
- known CAR-T targets are still high in the KO cells: CD276/B7-H3 175 TPM, EPHA2 556, CD70 364, CD44 881, PDPN 44, IL13RA2 50 (up), EGFR 12 (up), B4GALNT1 (GD2 synthase) 48. ROBO1 is the only one that drops (results/known_targets_check.csv, figures/known_targets). probably more useful for picking a co-target than the upregulated list

## secreted (TME)
- 15 secreted DE genes, 1 up (ADAMTS6) and 14 down
- down: CHI3L1 (YKL-40), POSTN, SPON1, SERPINE2, IGFBP5, THBS2, COL21A1, ADAMTS14, HSPG2, APOD, AZGP1 etc, mostly ECM
- KO cells might secrete less matrix and fewer myeloid recruiting factors (CHI3L1, POSTN)

## caveats / to check
- ROBO2 down with ROBO1. co-regulation? large deletion or regulatory element hit at the cut site? gRNA off-target? do genomic PCR + off-target check
- IFN genes can change from CRISPR/lentivirus handling or sorting, not just ROBO1. both groups went through the same process so probably ok
- 4 of the 12 up genes are unnamed/readthrough loci <7 TPM, could be quant artefacts
- XIST a bit lower in KO (log2FC -0.34), probably X inactivation drift in culture
- all replicates are from one line (BT241), can't generalize to GBM from this
- GENCODE v50 transcripts.fa had alt haplotype duplicates, removed before indexing
- trimmed fastqs deleted after quant (disk full), 01_qc.sh regenerates them

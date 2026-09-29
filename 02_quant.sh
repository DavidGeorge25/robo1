#!/bin/bash
# salmon index + quant against gencode v50 transcripts
# no genome decoy: only 24 gb ram on this mac

export MAMBA_ROOT_PREFIX=~/micromamba
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd ~/robo1_analysis
mkdir -p ref salmon

# reference files
gc=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50
[ -f ref/gencode.v50.transcripts.fa.gz ] || curl -s -o ref/gencode.v50.transcripts.fa.gz $gc/gencode.v50.transcripts.fa.gz
[ -f ref/gencode.v50.annotation.gtf.gz ] || curl -s -o ref/gencode.v50.annotation.gtf.gz $gc/gencode.v50.annotation.gtf.gz

# transcript -> gene table (tx, gene id, symbol, gene type, chr)
gzip -dc ref/gencode.v50.annotation.gtf.gz | perl -ne 'next unless /\ttranscript\t/; @f = split /\t/;
    ($g) = /gene_id "([^"]+)"/; ($t) = /transcript_id "([^"]+)"/; ($n) = /gene_name "([^"]+)"/; ($y) = /gene_type "([^"]+)"/;
    print "$t\t$g\t$n\t$y\t$f[0]\n"' > ref/tx2gene.tsv

# v50 transcripts.fa also has ~26k alt haplotype copies (e.g. 6 RPS18 genes on the MHC alts)
# keep only transcripts in the main gtf so reads aren't split between copies
cut -f1 ref/tx2gene.tsv > ref/keep.txt
gzip -dc ref/gencode.v50.transcripts.fa.gz | awk -F'|' 'NR == FNR {keep[$1]; next} /^>/ {p = (substr($1, 2) in keep)} p' ref/keep.txt - | gzip > ref/transcripts_primary.fa.gz
rm ref/keep.txt

# index, k=25 since reads are only 51 bp (default 31 is meant for 75 bp+)
[ -d ref/salmon_idx ] || salmon index -t ref/transcripts_primary.fa.gz -i ref/salmon_idx -k 25 --gencode -p 10

# quant, auto library type
for s in SS1-KO-1 SS2-KO-2 SS3-KO-3 SS4-AAVS1-1 SS5-AAVS1-2 SS6-AAVS1-3; do
    salmon quant -i ref/salmon_idx -l A \
        -1 trimmed/${s}_R1.fastq.gz -2 trimmed/${s}_R2.fastq.gz \
        --validateMappings --gcBias -p 10 -o salmon/$s 2> salmon/$s.log
done

# trimmed reads are ~7 gb and the disk is nearly full, remove once quantified (01_qc.sh remakes them)
rm -r trimmed

# multiqc again with salmon mapping rates added
multiqc qc/fastp salmon -o qc/multiqc -n multiqc_all --no-data-dir -f

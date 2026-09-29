#!/bin/bash
# trim + qc with fastp, then multiqc

export MAMBA_ROOT_PREFIX=~/micromamba
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

raw=~/robo1
cd ~/robo1_analysis
mkdir -p trimmed qc/fastp

samples="SS1-KO-1 SS2-KO-2 SS3-KO-3 SS4-AAVS1-1 SS5-AAVS1-2 SS6-AAVS1-3"

# fastp per sample
# nextseq 2000 (2-colour) so force polyG trimming, drop reads < 25 bp (salmon k)
for s in $samples; do
    d=${s%%-*}
    fastp -i $raw/$d/${s}_*_R1_001.fastq.gz -I $raw/$d/${s}_*_R2_001.fastq.gz \
        -o trimmed/${s}_R1.fastq.gz -O trimmed/${s}_R2.fastq.gz \
        --detect_adapter_for_pe --trim_poly_g --length_required 25 \
        -w 8 -h qc/fastp/$s.html -j qc/fastp/$s.json 2> qc/fastp/$s.log
done

# combined report
multiqc qc/fastp -o qc/multiqc -n multiqc_fastp --no-data-dir -f

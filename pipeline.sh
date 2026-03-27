#!/bin/bash
# BINF6110 Assignment 3 - Shotgun Metagenomics Pipeline
# De Filippis et al. 2019 (SRP126540)
# 3 vegan: SRR8146990, SRR8146961, SRR8146989
# 3 omnivore: SRR8146956, SRR8146969, SRR8146975

mkdir -p ~/A3_metagenomics/{raw_reads,trimmed,fastqc_results,kraken_output,bracken_output,results}

cd ~/A3_metagenomics/raw_reads
for sample in SRR8146990 SRR8146961 SRR8146989 SRR8146969 SRR8146975 SRR8146956; do
    fastq-dump ${sample} --split-files --gzip
done

for sample in SRR8146990 SRR8146961 SRR8146989 SRR8146969 SRR8146975 SRR8146956; do
    fastp -i ${sample}_1.fastq.gz -I ${sample}_2.fastq.gz \
        -o ~/A3_metagenomics/trimmed/${sample}_1_trim.fastq.gz \
        -O ~/A3_metagenomics/trimmed/${sample}_2_trim.fastq.gz \
        -q 20 -w 4 \
        -h ~/A3_metagenomics/fastqc_results/${sample}_fastp.html \
        -j ~/A3_metagenomics/fastqc_results/${sample}_fastp.json
    fastqc ~/A3_metagenomics/trimmed/${sample}_1_trim.fastq.gz \
           ~/A3_metagenomics/trimmed/${sample}_2_trim.fastq.gz \
           -o ~/A3_metagenomics/fastqc_results/ -t 4
    rm ${sample}_1.fastq.gz ${sample}_2.fastq.gz
    kraken2 --db ~/A3_metagenomics/database \
        --confidence 0.15 --threads 4 \
        --output /dev/null \
        --report ~/A3_metagenomics/kraken_output/${sample}.report \
        --paired \
        ~/A3_metagenomics/trimmed/${sample}_1_trim.fastq.gz \
        ~/A3_metagenomics/trimmed/${sample}_2_trim.fastq.gz
    rm ~/A3_metagenomics/trimmed/${sample}_1_trim.fastq.gz \
       ~/A3_metagenomics/trimmed/${sample}_2_trim.fastq.gz
done

cd ~/A3_metagenomics/bracken_output
for sample in SRR8146990 SRR8146961 SRR8146989 SRR8146969 SRR8146975 SRR8146956; do
    bracken -d ~/A3_metagenomics/database \
        -i ~/A3_metagenomics/kraken_output/${sample}.report \
        -o ${sample}.bracken \
        -w ${sample}_bracken.report \
        -r 300 -l S
done

kraken-biom SRR8146990_bracken.report SRR8146961_bracken.report \
    SRR8146989_bracken.report SRR8146969_bracken.report \
    SRR8146975_bracken.report SRR8146956_bracken.report \
    -o ~/A3_metagenomics/results/table.biom

multiqc ~/A3_metagenomics/fastqc_results/ -o ~/A3_metagenomics/fastqc_results/multiqc_report

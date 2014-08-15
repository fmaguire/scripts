#!/bin/bash
for id in Sample*;
    do cd $id;
    mkdir -p trimmomatic_output;
    cd trimmomatic_output;
    R1=../raw_illumina_reads/*R1_001.fastq;
    R2=../raw_illumina_reads/*R2_001.fastq;
    trimmomatic_cmd="java -jar /sequencing_projects/tools/Trimmomatic-0.32/trimmomatic-0.32.jar  \
    PE \
    -threads 14 \
    -trimlog ${id}_trimmomatic_log \
    $R1 \
    $R2 \
    ${id}_R1_paired_output \
    ${id}_R2_paired_output \
    ${id}_R1_unpaired_output \
    ${id}_R2_unpaired_output \
    ILLUMINACLIP:/sequencing_projects/tools/Trimmomatic-0.32/adapters/exeter_sequencing_adaptors.fasta:2:40:15 \
    LEADING:5 \
    TRAILING:5 \
    MINLEN:36";
    eval $trimmomatic_cmd
    cd ../..;
done;

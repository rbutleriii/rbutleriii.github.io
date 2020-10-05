---
layout: post
title: "Big batching your RNAseq"
author: "RRBIII"
categories: rnaseq
tags: [rnaseq,salmon,STAR]
image: pca-1.png
---


#### So did you do read alignment or pseudoalignment? Why not both?

Everyone will have an opinion on whether psuedoalignment is a good thing, so if you are going to use it, at least one collaborator will inevitably suggest you go back and try it the "other" way. Head them off and do both, since pseudoalignment costs very little compute.  

```sh
#!/usr/bin/env bash
exec 1>> command.log 2>&1
set -ex

##################################
# for running all of the trimming, alignment and alignment-free pipelines
# does not make use of unpaired outputs
# this script is set for a server with 100 threads and 256GB memory
##################################
# directory where everything will go
BASE_DIR=/path/to/base_dir
# directory of reads must have *_R1_001.fastq.gz and *_R2_001.fastq.gz filenames
READ_DIR=${BASE_DIR}/read_folder
# prior to running, the appropriate reference files must be generated in refs:
# an assembly gtf (must match version that built indexes below):
GTF=${BASE_DIR}/refs/Rattus_norvegicus.Rnor_6.0.96.gtf
# a salmon index:
SALMON_INDEX=${BASE_DIR}/refs/Rnor6.quasi_index_21
# a STAR index:
STAR_INDEX=${BASE_DIR}/refs/Rnor6_star
# there is a script in refs, modify for correct species and assembly
```

#### to trim or not to trim

With better and better sequence quality, this is making less of a difference, but if for instance you have custom adapters, probably a good idea. trimmomatic is also good, but why deal with java when fastp is pretty fast.

```sh
##################################
# trim reads
mkdir ${BASE_DIR}/trimmed
cd ${BASE_DIR}/trimmed
# use paired end reads from read_dir
# trying fastp, uses 2 worker threads per run (approx 250%-400% cpu)
# for quality trimming, use -5 and -3 with default window(4) and mean q(20)
trim(){
  out="$(basename -s '_R1_001.fastq.gz' $1)"
  R2="$(echo $1 | sed s/_R1_001.fastq.gz/_R2_001.fastq.gz/)"
  fastp \
    --thread 2\
    -i $1 \
    -I $R2 \
    -o ${out}_R1_trim.fastq.gz \
    -O ${out}_R2_trim.fastq.gz \
    --adapter_sequence AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
    --adapter_sequence_r2 GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG \
    -3 \
    -5 \
    -l 20 \
    -h ${out}.html \
    -j ${out}.json
}

export -f trim
parallel -j50 trim {} ::: ${READ_DIR}/*R1_001.fastq.gz
```

#### salmon

```sh
#################################
# Salmon pipeline
mkdir ${BASE_DIR}/salmon
cd ${BASE_DIR}/salmon
# for paired end data
# salmon counts quasi-mapping task
salmon_map(){
  out="$(basename -s '_R1_trim.fastq.gz' $1)"
  R2="$(echo $1 | sed 's/_R1_trim/_R2_trim/')"
  salmon --no-version-check quant \
    --seqBias \
    --gcBias \
    --validateMappings \
    -i $2 \
    -l A \
    -p 20 \
    -1 <(pigz -dc "$1") \
    -2 <(pigz -dc "$R2") \
    -o ${out}
}

# run loop
export -f salmon_map
parallel -j5 salmon_map {} $SALMON_INDEX \
  ::: ${BASE_DIR}/trimmed/*_R1_trim.fastq.gz
```

#### STAR using shared memory

An advantage for using STAR is that it allows you to put the reference into shared memory and run several samples without eating up a lot of memory. Note the `genomeLoad` options

```sh
#################################
# STAR-featureCounts pipeline
mkdir ${BASE_DIR}/bam
cd ${BASE_DIR}/bam
# for paired end data
# STAR mapping task
star_map(){
  out="$(basename -s '_R1_trim.fastq.gz' $1)"
  R2="$(echo $1 | sed 's/_R1_trim/_R2_trim/')"
  STAR \
    --runThreadN 20 \
    --genomeDir $2 \
    --genomeLoad LoadAndKeep \
    --readFilesCommand zcat \
    --outFileNamePrefix ${out}_star \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMmapqUnique 60 \
    --outSAMattrRGline ID:$out LB:lib1 PL:illumina PU:unit1 SM:$out \
    --outBAMsortingThreadN 10 \
    --outBAMsortingBinsN 10 \
    --limitBAMsortRAM 10000000000 \
    --outBAMcompression 10 \
    --bamRemoveDuplicatesType UniqueIdentical \
    --outFilterIntronMotifs RemoveNoncanonical \
    --outFilterMismatchNmax 2 \
    --outFilterScoreMinOverLread 0.30 \
    --outFilterMatchNminOverLread 0.30 \
    --alignSoftClipAtReferenceEnds No\
    --readFilesIn $1 $R2
}

# add to shared memory
STAR \
  --runThreadN 100 \
  --genomeDir $STAR_INDEX \
  --genomeLoad LoadAndExit

# run loop
export -f star_map
parallel -j5 star_map {} $STAR_INDEX ::: ${BASE_DIR}/trimmed/*_R1_trim.fastq.gz

# remove from shared memory
STAR \
  --genomeDir $STAR_INDEX \
  --genomeLoad Remove

# featureCounts for all in folder (max threads 64)
FC_OUTFILE="featureCounts"
featureCounts \
  -T 64 \
  -p \
  -B \
  -C \
  -D 2000 \
  -t exon \
  -g gene_id \
  -a $GTF \
  -o ${FC_OUTFILE}_counts.txt *.bam

########################################
# go home
cd ${BASE_DIR}

```

#### Next time
There it is, that script with 100 threads and 256GB memory can run a hundred samples easily overnight. In the future I can speed test the actual throughput.  



#!/usr/bin/env bash

# create environment
# mamba create -n ont "python<3.9" "samtools>1.10" raven-assembler minimap2 mappy pigz seqkit yacrd racon pip htslib bcftools pysam unimap
# conda activate ont
# python -m pip install "medaka==1.2.4"

# few simple helper functions
log() {
    NOW=$(date +"%m-%d %T")
    echo "[${NOW}] $@"
}

run() {
    echo "  CMD: $@"
    eval $@
}

# script to assemble ONT data
if [ -z "$2" ]; then
    echo "Usage: assemble.sh input-ont-reads.fastq.gz output-basename"
    exit 1
fi

# set some variables you will use later
FASTQ=$1
COV=100
MINLEN=5000
THREADS=4
SIZE=30400000
TMPDIR=$(uuidgen)
ASSEMBLY="$2.raven.racon.medaka.fasta"
SCRUBB="$2.reads-scrubbed.fastq.gz"

# make working dir
mkdir -p ${TMPDIR}

# filter for minimum length
log "Filtering ONT data by minimum read length ${MINLEN}"
cmd="seqkit seq --quiet -g --min-len ${MINLEN} ${FASTQ} | pigz -p ${THREADS} > ${TMPDIR}/len.fastq.gz"
run $cmd

if [ ! -f "$TMPDIR/len.fastq.gz" ]; then
    log "Read length filtering has failed"
    exit 1
fi

# calc coverage and subsample to COV
target_cov=$((COV*SIZE))
total=$(seqkit stats ${TMPDIR}/len.fastq.gz | grep -v '^file' | awk '{ print $5 }' | sed 's/,//g')
log "Filtered data is ${total} bp"
if [ "$total" -gt "$target_cov" ]; then
    propcov=$(bc <<< "scale=3; ${SIZE}/${total}"*$COV)
    log "Sampling data by ${propcov} to get to ${COV}X coverage"
    cmd="seqkit sample --quiet -s 12345 -p ${propcov} ${TMPDIR}/len.fastq.gz | pigz -p ${THREADS} > ${TMPDIR}/len.sampled.fastq.gz"
    run $cmd
else
    log "Data is less than ${COV}X, using all reads"
    ln -s ${TMPDIR}/len.fastq.gz ${TMPDIR}/len.sampled.fastq.gz
fi

if [ ! -f "$TMPDIR/len.sampled.fastq.gz" ]; then
    log "Sampling has failed, exiting"
    exit 1
fi

# now run through YACRD
log "Removing chimeric reads with YACRD"
cmd="minimap2 -x ava-ont -g 500 ${TMPDIR}/len.sampled.fastq.gz ${TMPDIR}/len.sampled.fastq.gz 2> /dev/null 1> ${TMPDIR}/overlaps.paf"
run $cmd
cmd="yacrd -i ${TMPDIR}/overlaps.paf -o ${TMPDIR}/report.yacrd -c 4 -n 0.4 scrubb -i ${TMPDIR}/len.sampled.fastq.gz -o ${TMPDIR}/len.sampled.scrubb.fastq.gz"
run $cmd

if [ ! -f "$TMPDIR/len.sampled.scrubb.fastq.gz" ]; then
    log "Chimera filtering with YACRD failed, exiting"
    exit 1
fi

# now assemble with Raven
log "Assembling chimera filtered reads with Raven"
cmd="raven -p 0 -t ${THREADS} ${TMPDIR}/len.sampled.scrubb.fastq.gz 2> /dev/null 1> ${TMPDIR}/raven.fasta"
run $cmd

if [ -f "raven.cereal" ]; then
    rm raven.cereal
fi

if [ ! -f "$TMPDIR/raven.fasta" ]; then
    log "de novo assembly with Raven failed, exiting"
    exit 1
fi
# should check for low coverage here.....

# now run racon
log "Error correcting assembly with Racon"
cmd="minimap2 -x map-ont -t ${THREADS} ${TMPDIR}/raven.fasta ${TMPDIR}/len.sampled.scrubb.fastq.gz 2> /dev/null 1> ${TMPDIR}/racon.paf"
run $cmd
cmd="racon -m 8 -x -6 -g -8 -w 500 -u -t ${THREADS} ${TMPDIR}/len.sampled.scrubb.fastq.gz ${TMPDIR}/racon.paf ${TMPDIR}/raven.fasta 2> /dev/null 1> ${TMPDIR}/raven.racon.fasta"
run $cmd

if [ ! -f "$TMPDIR/raven.racon.fasta" ]; then
    log "Racon error correcting failed, exiting"
    exit 1
fi

# now run medaka
log "Error correcting step 2 with Medaka"
cmd="medaka_consensus -f -i ${TMPDIR}/len.sampled.scrubb.fastq.gz -d ${TMPDIR}/raven.racon.fasta -o ${TMPDIR}/medaka -t ${THREADS} > ${TMPDIR}/medaka.log 2>&1"
run $cmd

if [ ! -f "${TMPDIR}/medaka/consensus.fasta" ]; then
    log "Medaka error correcting failed, exiting"
    exit 1
fi

# get final output and clean up
mv ${TMPDIR}/medaka/consensus.fasta ${ASSEMBLY}
mv ${TMPDIR}/len.sampled.scrubb.fastq.gz ${SCRUBB}
log "Evaluating some stats of assembly"
cmd="seqkit stats ${ASSEMBLY}"
run $cmd

log "Assembly written to: ${ASSEMBLY}"
log "ONT scrubbed reads: ${SCRUBB}"
log "Finished."
rm -r ${TMPDIR}
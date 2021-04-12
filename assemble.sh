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

# set some default variables you will use later
COV=100
MINLEN=5000
THREADS=4
SIZE=30400000
TMPDIR=$(uuidgen)

# parse command line options
while (( "$#" )); do
  case "$1" in
    -f|--fastq)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        FASTQ=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    -o|--output)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        OUT=$2
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    -c|--coverage)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        COV=$2
        shift 2
      fi
      ;;
    -m|--minlen)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        MINLEN=$2
        shift 2
      fi
      ;;
    -t|--threads)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        THREADS=$2
        shift 2
      fi
      ;;
    -g|--genomesize)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        SIZE=$2
        shift 2
      fi
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
  esac
done

# script to assemble ONT data
if [ -z "$FASTQ" ]; then
    echo "Usage: assemble.sh -f input-ont-reads.fastq.gz -o output-basename"
    echo " [optional arguments]"
    echo "    -t, --threads     Number of threads to use. Default: 4"
    echo "    -m, --minlen      Minumum read length. Default: 5000"
    echo "    -g, --genomesize  Estimated genome size. Default: 30400000"
    echo "    -c, --coverage    Coverage to use for assembly. Default: 100"
    exit 1
fi

if [ -z "$OUT" ]; then
    echo "Usage: assemble.sh -f input-ont-reads.fastq.gz -o output-basename"
    echo " [optional arguments]"
    echo "    -t, --threads     Number of threads to use. Default: 4"
    echo "    -m, --minlen      Minumum read length. Default: 5000"
    echo "    -g, --genomesize  Estimated genome size. Default: 30400000"
    echo "    -c, --coverage    Coverage to use for assembly. Default: 100"
    exit 1
fi

# tell user what options are
log "Running lr-asm of ${FASTQ} with -c ${COV} -m ${MINLEN} -g ${SIZE} -t ${THREADS}"

# set the output files
ASSEMBLY="${OUT}.raven.racon.medaka.fasta"
SCRUBB="${OUT}.reads-scrubbed.fastq.gz"

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
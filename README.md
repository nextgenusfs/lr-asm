# lr-asm
Long read genome assembly scripts

#### Setup environment
1. Install miniconda3 (mac os link below)
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
echo 'alias miniconda="eval $(${HOME}/miniconda/bin/conda shell.bash hook)"' >> ~/.bash_profile
source ~/.bash_profile
```
2. Activate miniconda and install new environment
```
$ miniconda
$ conda create -n ont "python<3.9" "samtools>1.10" raven-assembler minimap2 mappy pigz seqkit yacrd racon pip htslib bcftools pysam
```
3. Activate environment and install medaka
```
$ conda activate ont
$ python -m pip install "medaka==1.2.4"
```
#### Running the assembly script
```
$ assemble.sh -f ont-reads.fastq.gz -o basename
```


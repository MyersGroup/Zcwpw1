---
title: "Zcwpw1 Project"
output: html_document
---




```{r setup, include=FALSE}

# Load required software dependencies

export PATH=/homes/wells/saxony/single-cell/sequencing/software/fastqc/FastQC/:$PATH
export PATH=/homes/wells/saxony/single-cell/sequencing/software/bedtools/bedtools2/bin/:$PATH
export PATH=/homes/wells/saxony/single-cell/sequencing/software/samtools/samtools-1.7/:$PATH
export LD_LIBRARY_PATH=/homes/wells/saxony/single-cell/sequencing/software/htslib_1.9/lib/:$LD_LIBRARY_PATH
export PATH=/homes/wells/saxony/single-cell/sequencing/software/:$PATH
export PATH=/homes/wells/saxony/single-cell/sequencing/software/MAPeakCaller/:$PATH
export PATH=/homes/wells/saxony/single-cell/sequencing/software/bwa-0.7.17/:$PATH
export PATH=/homes/wells/saxony/single-cell/sequencing/software/bowtie/bowtie-1.2.2-linux-x86_64/:$PATH
export PATH=/homes/wells/saxony/single-cell/sequencing/software/bin/:$PATH # bedops


```


Details of the original files, the ID to description mappings are included in Sample_Manifest.txt

```{r}

# The directory containing fastq files for each group and the genome specification is contained in the file config.yml

source ~/saxony/anaconda3/bin/activate
conda init

# Map reads to genome etc.
snakemake --cores 15 -npr --config GROUP="659233"
snakemake --cores 15 -npr --config GROUP="538916"
snakemake --cores 15 -npr --config GROUP="594404"
snakemake --cores 15 -npr --config GROUP="Altemose2015"
snakemake --cores 15 -npr --config GROUP="Dmc1_r1"

# Call Peaks
snakemake --cores 15 --snakefile Snakefile_peaks -npr

# make enrichment profile plots (& center/strand by motif if possible)
snakemake --cores 15 --snakefile snakemake_profileplot -npr

# force-call
snakemake --cores 15 --snakefile Snakefile_forcecall -npr

# plot profile over gene
metagene_plot.sh

```


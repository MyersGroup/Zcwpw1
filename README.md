# Zcwpw1 Project

Zcwpw1 is co-expressed with Prdm9 and has domains that bind to both H3K4me3 and H3K36me3 (which are desposited by Prdm9 and required for its role in recombination.)

Software dependencies are managed as a Conda environment:

```{bash}
source ~/saxony/anaconda3/bin/activate
conda init

conda env create -f environment.yml
source activate zcwpw1
```

Details of the original FASTQ files, and ID-to-description mappings are included in [Sample_Manifest.Rmd](analysis/Sample_Manifest.Rmd)

Most of the analyses are arranged as snakemake pipelines for reproducibility.

The directory containing fastq files for each group and the genome specification is contained in the file [config.yml](pipelines/config.yml)

```{bash}
# Map reads to genome etc.
snakemake --cores 15 -npr --config GROUP="659233"
snakemake --cores 15 -npr --config GROUP="538916"
snakemake --cores 15 -npr --config GROUP="594404"
snakemake --cores 15 -npr --config GROUP="Altemose2015"
snakemake --cores 15 -npr --config GROUP="Dmc1_r1"

# Call Peaks
snakemake --cores 15 --snakefile pipelines/Snakefile_peaks -npr

# make enrichment profile plots (& center/strand by motif if possible)
snakemake --cores 15 --snakefile pipelines/snakemake_profileplot -npr

# force-call
snakemake --cores 15 --snakefile pipelines/Snakefile_forcecall -npr

# plot profile over gene
metagene_plot.sh
```


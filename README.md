# Zcwpw1 Project

Zcwpw1 is co-expressed with Prdm9 and has domains that bind to both H3K4me3 and H3K36me3 (which are desposited by Prdm9 and required for its role in recombination.)

Software dependencies are managed as a Conda environment detailed in [setup.sh](setup.sh)

Details of the original FASTQ files, and ID-to-description mappings are included in [Sample_Manifest.Rmd](analysis/Sample_Manifest.Rmd)

The directory containing fastq files for each group and the genome specification is contained in the file [config.yml](pipelines/config.yml)

Most of the analyses are arranged as snakemake pipelines for reproducibility:

```{bash}
# Map reads to genome etc.
snakemake --snakefile pipelines/Map_Reads.py -npr --dag --forceall --config GROUP="Dmc1_r1" | dot -Tpdf > Map_Reads_dag.pdf
snakemake --cores 15 --snakefile pipelines/Map_Reads.py -npr --config GROUP="659233"
snakemake --cores 15 --snakefile pipelines/Map_Reads.py -npr --config GROUP="538916"
snakemake --cores 15 --snakefile pipelines/Map_Reads.py -npr --config GROUP="594404"
snakemake --cores 15 --snakefile pipelines/Map_Reads.py -npr --config GROUP="NA15"
snakemake --cores 15 --snakefile pipelines/Map_Reads.py -npr --config GROUP="Dmc1_r1"

# Call Peaks
snakemake --cores 15 --snakefile pipelines/Call_Peaks.py -npr

# make enrichment profile plots (& center/strand by motif if possible)
snakemake --cores 15 --snakefile pipelines/Plot_Profile2.py -npr

snakemake --cores 15 --snakefile pipelines/Plot_Heatmap.py -npr

# force-call
# specify sample pairs to force call in sample_pairings.py
# specify locations to call at in config.yml
snakemake --cores 15 --snakefile pipelines/Force_Call_Peaks.py -npr

snakemake --cores 15 --snakefile pipelines/wgs.py -npr

# create DMC1 profile plot
cd dmc1
snakemake --cores 15 --snakefile ../pipelines/dmc1.py -npr

#Create Allele Specificify Plot
Allele_specificity.sh

# plot profile over gene
metagene_plot.sh
```

Other analyses are Rnotebooks, and require the [preprocess.sh](preprocess.sh) script to be run to generate the required files.



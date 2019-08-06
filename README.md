# Zcwpw1 Project

Zcwpw1 is co-expressed with Prdm9 and has domains that bind to both H3K4me3 and H3K36me3 (which are desposited by Prdm9 and required for its role in recombination.)

Software dependencies are managed as a Conda environment:

```{bash}
#conda config --set auto_activate_base false

source ~/saxony/anaconda3/bin/activate
conda init

conda env create -f environment.yml
source activate zcwpw1

conda deactivate
conda remove --name zcwpw1 --all


conda install -c bioconda ucsc-bedgraphtobigwig

echo ".libPaths( c( '~/R/' , .libPaths() ) )" > ~/.Rprofile
install.packages("remotes")
remotes::install_github("myersgroup/MotifFinder")

#conda install openssl=1.0
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
chmod 777 wigToBigWig
./wigToBigWig

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
chmod 777 bigWigToBedGraph


# conda env remove --name zcwpw1
conda install -c conda-forge ripgrep


# install bwtool
git clone https://github.com/CRG-Barcelona/bwtool.git
git clone https://github.com/CRG-Barcelona/libbeato.git
cd libbeato
git checkout 0c30432
./configure --prefix=$HOME CFLAGS="-g -O0 -I${HOME}/include" LDFLAGS=-L${HOME}/lib
make
make install
cd ../bwtool/
./configure --prefix=$HOME CFLAGS="-g -O0 -I${HOME}/include" LDFLAGS=-L${HOME}/lib
make
make install

```

Details of the original FASTQ files, and ID-to-description mappings are included in [Sample_Manifest.Rmd](analysis/Sample_Manifest.Rmd)

Most of the analyses are arranged as snakemake pipelines for reproducibility.

The directory containing fastq files for each group and the genome specification is contained in the file [config.yml](pipelines/config.yml)

```{bash}
# Map reads to genome etc.
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
snakemake --cores 15 --snakefile pipelines/Force_Call_Peaks.py -npr

# plot profile over gene
metagene_plot.sh
```


Create Allele Specificify Plot
```{bash}
Rscript pipelines/MultiProfilePlot.R bwplots/AlleleSpecificity.pdf 'Zcwpw1 binds Prdm9 sites in an allele specific manner' 10 8 \
bwprofilesNorm/WTCHG_538916_223180_VS_WTCHG_538916_221156_AT_SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL.profile \
bwprofilesNorm/WTCHG_538916_223180_VS_WTCHG_538916_221156_AT_SingleBasePeaks.NA15-SRR5627145_AND_NA15-SRR5627144_vs_NA15-SRR5627143.p0.000001.sep250.ALL.profile \
bwprofilesNorm/WTCHG_538916_224192_VS_WTCHG_538916_221156_AT_SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL.profile \
bwprofilesNorm/WTCHG_538916_224192_VS_WTCHG_538916_221156_AT_SingleBasePeaks.NA15-SRR5627145_AND_NA15-SRR5627144_vs_NA15-SRR5627143.p0.000001.sep250.ALL.profile \
'1) Zcwpw1 cotransfected with human Prdm9, at human Prdm9 sites' \
'2) Zcwpw1 cotransfected with human Prdm9, at chimp Prdm9 sites' \
'3) Zcwpw1 cotransfected with chimp Prdm9, at human Prdm9 sites' \
'4) Zcwpw1 cotransfected with chimp Prdm9, at chimp Prdm9 sites'

# Use Motif Centered and Stranded for Human allele
Rscript pipelines/MultiProfilePlot.R bwplots/AlleleSpecificityMCT.pdf 'Zcwpw1 binds Prdm9 sites in an allele specific manner (Motif Centered and Stranded MCS)' 11 8 \
bwprofilesNorm/WTCHG_538916_223180_VS_WTCHG_538916_221156_AT_SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL_MotifCenteredStranded.profile \
bwprofilesNorm/WTCHG_538916_223180_VS_WTCHG_538916_221156_AT_SingleBasePeaks.NA15-SRR5627145_AND_NA15-SRR5627144_vs_NA15-SRR5627143.p0.000001.sep250.ALL.profile \
bwprofilesNorm/WTCHG_538916_224192_VS_WTCHG_538916_221156_AT_SingleBasePeaks.NA15-SRR5627146_AND_NA15-SRR5627147_vs_NA15-SRR5627143.p0.000001.sep250.ALL_MotifCenteredStranded.profile \
bwprofilesNorm/WTCHG_538916_224192_VS_WTCHG_538916_221156_AT_SingleBasePeaks.NA15-SRR5627145_AND_NA15-SRR5627144_vs_NA15-SRR5627143.p0.000001.sep250.ALL.profile \
'1) Zcwpw1 cotransfected with human Prdm9, at human Prdm9 MCS sites' \
'2) Zcwpw1 cotransfected with human Prdm9, at chimp Prdm9 sites' \
'3) Zcwpw1 cotransfected with chimp Prdm9, at human Prdm9 MCS sites' \
'4) Zcwpw1 cotransfected with chimp Prdm9, at chimp Prdm9 sites'
```

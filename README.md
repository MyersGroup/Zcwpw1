# Zcwpw1 Project

Zcwpw1 is co-expressed with Prdm9 and has domains that bind to both H3K4me3 and H3K36me3 (which are desposited by Prdm9 and required for its role in recombination.)

Software dependencies are managed as a Conda environment:

```{bash}
#conda config --set auto_activate_base false

source ~/saxony/anaconda3/bin/activate
conda init

conda env create -f environment.yml
source activate zcwpw1
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
snakemake --cores 15 --snakefile pipelines/Map_Reads.py -npr --config GROUP="Altemose2015"
snakemake --cores 15 --snakefile pipelines/Map_Reads.py -npr --config GROUP="Dmc1_r1"

# Call Peaks
snakemake --cores 15 --snakefile pipelines/Call_Peaks.py -npr

# make enrichment profile plots (& center/strand by motif if possible)
snakemake --cores 15 --snakefile pipelines/Plot_Profile2.py -npr

# force-call
snakemake --cores 15 --snakefile pipelines/Force_Call_Peaks.py -npr

# plot profile over gene
metagene_plot.sh
```

#conda config --set auto_activate_base false

source ~/saxony/anaconda3/bin/activate
conda init

conda env create -f environment.yml
source activate zcwpw1

conda deactivate
conda remove --name zcwpw1 --all

#conda install -c bioconda ucsc-liftover
#conda install -c bioconda crossmap
#pip3 install CrossMap --user "No module named bx"
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
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

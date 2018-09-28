#!/bin/bash

# unload potentially conflicting anaconda instances
{ # try
    module unload python-anaconda2 &&
    module unload python-anaconda3
} || { # catch
    echo 'module unloading failed: maybe module does not exist'
}


# install miniconda for local independent package management
wget http://repo.continuum.io/miniconda/Miniconda2-4.3.21-Linux-x86_64.sh -O miniconda.sh
mkdir dependencies
chmod 775 miniconda.sh
chmod 775 dependencies
bash miniconda.sh -b -p ./dependencies/miniconda
rm miniconda.sh

# activate conda virtual environment
source ./dependencies/miniconda/bin/activate

# add bioconda and r channel for easy dependency installations
conda config --add channels r
conda config --add channels bioconda
pip install screed
conda install -c bioconda mothur bioconductor-ggtree fasttree blast -y
conda install -c conda-forge readline -y

# install pyflow for automated task management
wget https://github.com/Illumina/pyflow/releases/download/v1.1.17/pyflow-1.1.17.tar.gz
pip install pyflow-1.1.17.tar.gz
rm pyflow-1.1.17.tar.gz

# install MeFiT and dependencies
conda install numpy click jellyfish HTSeq

wget http://best.snu.ac.kr/casper/program/casper_v0.8.2.tar.gz
tar -zxvf casper_v0.8.2.tar.gz
rm casper_v0.8.2.tar.gz
cd casper_v0.8.2 && make
mv casper ../dependencies/miniconda/bin/
cd .. && rm -r casper_v0.8.2

wget https://raw.githubusercontent.com/nisheth/MeFiT/master/mefit
chmod 772 mefit
mv mefit dependencies/miniconda/bin/

# install Silva databases
mkdir -p dependencies/silva
wget https://mothur.org/w/images/7/71/Silva.seed_v132.tgz
tar -zxvf Silva.seed_v132.tgz -C dependencies/silva/

wget https://mothur.org/w/images/3/32/Silva.nr_v132.tgz
tar -zxvf Silva.nr_v132.tgz -C dependencies/silva/
rm -r Silva*

# install TaxAss
git clone https://github.com/McMahonLab/TaxAss.git
unzip TaxAss/FreshTrain-files/FreshTrain18Aug2016.zip
mv TaxAss dependencies/
mv FreshTrain18Aug2016 dependencies/TaxAss/

# install vsearch
wget https://github.com/torognes/vsearch/archive/v2.8.5.tar.gz
tar zxvf v2.8.5.tar.gz
cd vsearch-2.8.5
./autogen.sh
./configure --prefix=$HOME
make
make install
cp bin/vsearch ../dependencies/miniconda/bin/
cd ..
rm -r vsearch*
rm v2.8.5.tar.gz

# Omics 16S
Base workflow for setting up OTU analysis from Illumina fastq(s).  It follows the [Mothur MiSeq SOP](https://mothur.org/wiki/MiSeq_SOP) except for the initial quality control steps.  High-quality read length and abundance have been demonstrated to be primary factors in avoiding spurious and/or inflated OTU classification, so special emphasis was placed on those steps.

  * [Quality-filtering vastly improves diversity estimates from Illumina amplicon sequencingi 2013](https://www.nature.com/nmeth/journal/v10/n1/full/nmeth.2276.html)
  * [Accuracy of microbial community diversity estimated by closed- and open-reference OTUs 2017](https://peerj.com/articles/3889/)


It leverages third party tools:

  * [MeFit 2016](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1358-1): a merging and filtering tool for Illumina paired-end reads, designed specifically for 16S rRNA amplicaon sequencing data
  * [Mothur 2009](http://aem.asm.org/content/75/23/7537.full): software for describing and comparing microbial communities

## Setup
1. To be able to build [CASPER](http://best.snu.ac.kr/casper/index.php?name=manual), you'll need to make sure the g++ compiler and boost libraries are installed

2. Build the local conda environment with dependencies
> $./build

3. Activate the root environment
> $source <path_to_project>/dependencies/miniconda/bin/activate

4. Run and build on scripts at will!
> $python runner.py <run_dp>

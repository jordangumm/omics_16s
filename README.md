# Omics 16S
Automated workflow for OTU classification using shotgun data from Illumina Miseq fastq(s).  It follows the [Mothur MiSeq SOP](https://mothur.org/wiki/MiSeq_SOP) except for the initial quality control steps.  High-quality read length and abundance have been demonstrated to be primary factors in avoiding spurious and/or inflated OTU classification, so special emphasis was placed on those steps.  It also attempts to refine freshwater OTU classifications using TaxAss.

  * [Quality-filtering vastly improves diversity estimates from Illumina amplicon sequencingi 2013](https://www.nature.com/nmeth/journal/v10/n1/full/nmeth.2276.html)
  * [Accuracy of microbial community diversity estimated by closed- and open-reference OTUs 2017](https://peerj.com/articles/3889/)


It leverages third party tools and databases:

  * [MeFit 2016](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1358-1): a merging and filtering tool for Illumina paired-end reads, designed specifically for 16S rRNA amplicaon sequencing data
  * [Mothur 2009](http://aem.asm.org/content/75/23/7537.full): software for describing and comparing microbial communities
  * [Silva reference files v128 2016](https://mothur.org/wiki/Silva_reference_files): 16S rRNA seed database and  sequence/taxonomy references
  * [TaxAss 2017](https://www.biorxiv.org/content/early/2017/11/05/214288): fine-scale taxonomic assignment for freshwater datasets (by default)

## Setup
1. To be able to build [CASPER](http://best.snu.ac.kr/casper/index.php?name=manual), you'll need to make sure the g++ compiler and boost libraries are installed

2. Build the local conda environment with dependencies
> $./build

3. Activate the root environment
> $source <path_to_project>/dependencies/miniconda/bin/activate


## Run Analysis
Use the runner script against a sequencing run directory of fastq(s)
> $python runner.py <run_dp>

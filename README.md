# Omics 16S
Base workflow for setting up OTU analysis from Illumina fastq(s)

## Setup

1. Build the local conda environment with dependencies
> $./build

2. Activate the root environment
> $source <path_to_project>/dependencies/miniconda/bin/activate

3. Run and build on scripts at will!
> $python preanalysis.py <run_dp>

4. Install new dependencies
> $conda install bwa

5. Save environment dependencies for later automated build
> $conda env export > environment.yml

6. Recreate saved environment
> $conda env create -f environment.yml

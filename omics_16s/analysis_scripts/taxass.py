import sys
reload(sys)
sys.setdefaultencoding('utf8')

import os
import random
import string

import click
from pyflow import WorkflowRunner
from subprocess import call

class TaxAss(WorkflowRunner):
    def __init__(self, fasta_fp, taxass_dp, silva_dp, num_cpu):
        self.fasta_fp = fasta_fp

        self.specific_fasta = os.path.join(taxass_dp, 'FreshTrain18Aug2016.fasta')
        self.specific_taxa = os.path.join(taxass_dp, 'FreshTrain18Aug2016.taxonomy')
        if not os.path.exists(self.fasta_ref) or not os.path.exists(self.taxonomy_ref):
            sys.exit('[ERROR]: Taxass file not found in {}'.format(taxass_dp))

        self.general_fasta = os.path.join(silva_dp, 'silva.nr_
        self.general_taxa =
        self.num_cpu=num_cpu


    def workflow(self):
        """ method invoked on class instance run call """
        pass


@click.command()
@click.argument('read_dp')
@click.argument('taxass_dp')
@click.argument('silva_dp')
@click.option('-num_cpu', type=click.INT, default=2)
def analysis(read_dp, taxass_dp, silva_dp):
    """ TaxAss Analysis Management """
    log_output_dp = os.path.join(analysis_dp, 'logs')

    preanalyzer = Analyzer(read_dp=read_dp,
                           taxass_dp=taxass_dp,
                           silva_dp=silva_dp,
                           num_cpu=num_cpu)
    preanalyzer.run(mode='local', dataDirRoot=log_output_dp)


if __name__ == "__main__":
    analysis()

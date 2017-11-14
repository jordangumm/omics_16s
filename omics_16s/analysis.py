import sys
reload(sys)
sys.setdefaultencoding('utf8')

import os
import random
import string

import click
from pyflow import WorkflowRunner
from subprocess import call

class Analyzer(WorkflowRunner):
    def __init__(self, analysis_dp, num_cpu):
        self.analysis_dp = analysis_dp
        self.read_dp = os.path.join(analysis_dp, 'reads')
        self.mothur_fasta_fp = os.path.join(self.read_dp, 'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta')
        if not os.path.exists(self.read_dp) or not os.path.exists(self.mothur_fasta_fp):
            sys.exit('[ERROR]: No or missing analysis reads, preanalysis must have failed!')

        self.mothur_dp = os.path.join(analysis_dp, 'mothur')
        self.num_cpu=num_cpu

        full_dp = os.path.dirname(os.path.abspath(__file__))
        self.otu_generator_fp = os.path.join(full_dp, 'analysis_scripts', 'otu_generator.py')
        self.taxass_fp = os.path.join(full_dp, 'analysis_scripts', 'taxass.py')

        self.dependencies_dp = full_dp.replace('omics_16s/omics_16s', 'omics_16s/dependencies')


    def workflow(self):
        """ method invoked on class instance run call """
        self.addTask("otu_generator", command=['python', self.otu_generator_fp, '--num_cpu', num_cpu,
                                                                  self.read_dp, self.dependencies_dp])
        self.addTask("taxass_analysis", command=['python', self.taxass_fp, '--num_cpu', num_cpu,
                                   self.mothur_fasta_fp, os.path.join(self.dependencies_dp, 'silva')],
                                        dependencies=['otu_generator',])


@click.command()
@click.argument('analysis_dp')
@click.option('--num_cpu', '-n', default=2)
def analysis(analysis_dp, num_cpu):
    """ Analysis Management """
    log_output_dp = os.path.join(analysis_dp, 'logs')

    preanalyzer = Analyzer(analysis_dp=analysis_dp, num_cpu=num_cpu)
    preanalyzer.run(mode='local', dataDirRoot=log_output_dp)


if __name__ == "__main__":
    analysis()

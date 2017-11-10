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
    def __init__(self, analysis_dp, num_cpu=1):
        self.analysis_dp = analysis_dp
        self.read_dp = os.path.join(analysis_dp, 'reads')
        if not os.path.exists(self.read_dp):
            sys.exit('[ERROR]: No analysis reads, preanalysis must have failed!')
        self.num_cpu=num_cpu


    def workflow(self):
        """ method invoked on class instance run call """
        pass


@click.command()
@click.argument('analysis_dp')
def analysis(analysis_dp):
    """ Pre-Analysis Management

    Sets up Pyflow WorkflowRunner with randomized log output string to avoid 
    possible output collision if ran multiple times.

    Arguments:
    run_dp -- String path to run directory to pre-analyze
    """
    log_output_dp = os.path.join(analysis_dp, 'logs')

    preanalyzer = Analyzer(analysis_dp=analysis_dp)
    preanalyzer.run(mode='local', dataDirRoot=log_output_dp)


if __name__ == "__main__":
    analysis()

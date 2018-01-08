import sys
reload(sys)
sys.setdefaultencoding('utf8')

import os
import random
import string

import click
from pyflow import WorkflowRunner
from subprocess import call

class FastTreeRunner(WorkflowRunner):
    def __init__(self, fasta_fp, analysis_dp, num_cpu):
        self.output_dp = os.path.join(analysis_dp, 'fasttree')
        if not os.path.exists(self.output_dp):
            os.makedirs(self.output_dp)

        self.fasta_fp = fasta_fp
        #os.path.join(analysis_dp, 'mothur',
        #     'stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta')
        if not os.path.exists(self.fasta_fp):
            sys.exit('mothur files not available for taxass analysis')


    def workflow(self):
        """ method invoked on class instance run call """
        cmd = "FastTree -nt < {} > {}".format(self.fasta_fp, os.path.join(self.output_dp, 'jukes-cantor+cat.tre'))
        self.addTask('jukes', command=cmd)

        cmd = "FastTree -gtr -nt < {} > {}".format(self.fasta_fp, os.path.join(self.output_dp, 'gtr+cat.tre'))
        self.addTask('gtr', command=cmd, dependencies=['jukes',])





@click.command()
@click.argument('fasta_fp')
@click.argument('analysis_dp')
@click.option('--num_cpu', '-n', type=click.INT, default=4)
def analysis(fasta_fp, analysis_dp, num_cpu):
    """ Runs mothur fasta output against FastTree workflow

    \b
    Arguments:
        fasta_fp:    file path to fasta file to run through fasttree
        analysis_dp: directory path to active analysis directory
    """
    log_output_dp = os.path.join(analysis_dp, 'logs', 'fasttree')

    runner = FastTreeRunner(fasta_fp=fasta_fp,
                            analysis_dp=analysis_dp,
                            num_cpu=num_cpu)
    runner.run(mode='local', dataDirRoot=log_output_dp)

    return None



if __name__ == "__main__":
    analysis()

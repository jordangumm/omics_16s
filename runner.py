import sys
import os
import random
import string

import click
from pyflow import WorkflowRunner
from subprocess import call

class Runner(WorkflowRunner):
    def __init__(self, run_dp):
        self.run_dp = run_dp

    def workflow(self):
        """ method invoked on class instance run call """
        full_dp = os.path.dirname(os.path.abspath(__file__))
        self.addTask("preanalysis", command=['python', os.path.join(full_dp, 'omics_16s/preanalysis.py'), self.run_dp])
        self.addTask("analysis", command=['python', os.path.join(full_dp, 'omics_16s/analysis.py'),
                                                    os.path.join(self.run_dp, 'bioinfo', 'analysis')],
                                 dependencies=['preanalysis',])


@click.command()
@click.argument('run_dp')
@click.option('--flux/--no-flux', default=False)
@click.option('--account', '-a', default='lsa_fluxm')
@click.option('--ppn', '-p', default=8)
@click.option('--mem', '-m', default='20gb')
@click.option('--walltime', '-w', default='2:00:00')
def pre_analysis(run_dp, flux, account, ppn, mem, walltime):
    """ Pre-Analysis Management

    Sets up Pyflow WorkflowRunner and launches locally by default or via flux

    Arguments:
    run_dp -- String path to run directory to pre-analyze
    """
    log_output_dp = os.path.join(run_dp, 'bioinfo', 'runner')

    if flux:
        call('echo "python runner.py {}" | qsub -N omics_16s -A {} -q fluxm -l nodes=1:ppn={},mem={},walltime={}'.format(
                                                                              run_dp, account, ppn, mem, walltime), shell=True)
    else:
        runner = Runner(run_dp=run_dp)
        runner.run(mode='local', dataDirRoot=log_output_dp)


if __name__ == "__main__":
    pre_analysis()
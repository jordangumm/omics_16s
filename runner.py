import sys
import os
import random
import string

import click
from pyflow import WorkflowRunner
from subprocess import call

class Runner(WorkflowRunner):
    def __init__(self, run_dp, num_cpu):
        full_dp = os.path.dirname(os.path.abspath(__file__))
        self.preanalysis_fp = os.path.join(full_dp, 'omics_16s', 'preanalysis.py')
        self.analysis_fp = os.path.join(full_dp, 'omics_16s', 'analysis.py')
 
        self.silva_dp = os.path.join(full_dp, 'dependencies', 'silva')
        self.taxass_dp = os.path.join(full_dp, 'dependencies', 'TaxAss')

        self.run_dp = run_dp
        self.analysis_dp = os.path.join(run_dp, 'bioinfo', 'analysis')
        self.num_cpu = num_cpu

    def workflow(self):
        """ method invoked on class instance run call """
        self.addTask("preanalysis", command=['python', self.preanalysis_fp, self.run_dp, self.analysis_dp])
        self.addTask("analysis", command=['python', self.analysis_fp, '--num_cpu', self.num_cpu, self.analysis_dp], dependencies=['preanalysis',])


@click.command()
@click.argument('run_dp')
@click.option('--flux/--no-flux', default=False)
@click.option('--queue', '-q', default='fluxod')
@click.option('--account', '-a')
@click.option('--ppn', '-p', default=4)
@click.option('--mem', '-m', default='20gb')
@click.option('--walltime', '-w', default='2:00:00')
def runner(run_dp, flux, queue, account, ppn, mem, walltime):
    """ Analysis Workflow Management

    Sets up Pyflow WorkflowRunner and launches locally by default or via flux

    Arguments:
    run_dp -- String path to run directory to pre-analyze
    """
    log_output_dp = os.path.join(run_dp, 'bioinfo', 'logs', 'runner')

    if flux:
        full_dp = os.path.dirname(os.path.abspath(__file__))
        activate = 'source {}'.format(os.path.join(full_dp, 'dependencies', 'miniconda', 'bin', 'activate'))
        runner_fp = os.path.join(full_dp, 'runner.py')
        if queue == 'fluxod':
            qsub = 'qsub -N omics_16s -A {} -q {} -l nodes=1:ppn={}:largemem,mem={},walltime={}'.format(
                                                             account, queue, ppn, mem, walltime)
        else:
            qsub = 'qsub -N omics_16s -A {} -q {} -l nodes=1:ppn={},mem={},walltime={}'.format(
                                                             account, queue, ppn, mem, walltime)
        call('echo "{} && python {} {}" | {}'.format(activate, runner_fp, run_dp, qsub), shell=True)
    else:
        workflow_runner = Runner(run_dp=run_dp, num_cpu=ppn)
        workflow_runner.run(mode='local', dataDirRoot=log_output_dp)


if __name__ == "__main__":
    runner()

import sys
import os
import random
import string

import click
from pyflow import WorkflowRunner

class Runner(WorkflowRunner):
    def __init__(self, run_dp):
        self.run_dp = run_dp

    def workflow(self):
        """ method invoked on class instance run call """
        self.addTask("preanalysis", command=['python', 'preanalysis.py', self.run_dp])
        self.addTask("analysis", command=['python', 'analysis.py', os.path.join(self.run_dp, 'bioinfo', 'analysis')],
                                 dependencies=['preanalysis',])


@click.command()
@click.argument('run_dp')
def pre_analysis(run_dp):
    """ Pre-Analysis Management

    Sets up Pyflow WorkflowRunner with randomized log output string to avoid 
    possible output collision if ran multiple times.

    Arguments:
    run_dp -- String path to run directory to pre-analyze
    """
    log_output_dp = os.path.join(run_dp, 'bioinfo', 'runner')

    runner = Runner(run_dp=run_dp)
    runner.run(mode='local', dataDirRoot=log_output_dp)


if __name__ == "__main__":
    pre_analysis()

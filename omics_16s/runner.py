import sys
import os
import random
import string

import click
from pyflow import WorkflowRunner

class PreAnalyzer(WorkflowRunner):
    def __init__(self, run_dp):
        self.run_dp = run_dp

    def workflow(self):
        """ method invoked on class instance run call """
        self.addTask("flowcellqc",
                     command=['echo', '"test flowcellqc output"'])


@click.command()
@click.argument('run_dp')
def pre_analysis(run_dp):
    """ Pre-Analysis Management

    Sets up Pyflow WorkflowRunner with randomized log output string to avoid 
    possible output collision if ran multiple times.

    Arguments:
    run_dp -- String path to run directory to pre-analyze
    """
    randstr = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output_dp = os.path.join(run_dp, "preanalysis{}".format(randstr))

    preanalyzer = PreAnalyzer(run_dp=run_dp)
    preanalyzer.run(mode='local', dataDirRoot=log_output_dp)


if __name__ == "__main__":
    pre_analysis()

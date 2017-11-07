import sys
import os
import random
import string

import click
from pyflow import WorkflowRunner

class PreAnalyzer(WorkflowRunner):
    def __init__(self, run_dp):
        self.run_dp = run_dp
        self.samples = self.get_sample_dict()

    def get_sample_dict(self):
        sids = {}
        for fastq in os.listdir(self.run_dp):
            fastq_fp = os.path.join(self.run_dp, fastq.replace('.gz',''))
            if 'fastq' not in fastq: continue
            sid = '_'.join([x for x in fastq.split('_')][:-3])

            if sid not in sids.keys():
                if '_R1_' in fastq: sids[sid] = {'r1': fastq_fp}
                elif '_R2_' in fastq: sids[sid] = {'r2': fastq_fp}
                else: sys.exit('[ERROR]: {} not an R1 or R2 fastq file!'.format(fastq_fp))
            else:
                if '_R1_' in fastq: sids[sid]['r1'] = fastq_fp
                elif '_R2_' in fastq: sids[sid]['r2'] = fastq_fp
                else: sys.exit('[ERROR]: {} not an R1 or R2 fastq file!'.format(fastq_fp))
        return sids

    def get_gunzip_cmd(self):
        gunzip_cmd = ['gunzip']
        for fastq in os.listdir(self.run_dp):
            if 'fastq' not in fastq: continue
            if '.gz' in fastq:
                gunzip_cmd.append(os.path.join(self.run_dp, fastq))
        if len(gunzip_cmd) == 1: return ['echo', 'no files to ungzip']
        return gunzip_cmd

    def workflow(self):
        """ method invoked on class instance run call """
        self.addTask("gunzip", command=self.get_gunzip_cmd())
        for i, sid in enumerate(self.samples.keys()):
            self.addTask("mefit_{}".format(i),
                          command=['mefit', '-s', os.path.join('tmp', sid),
                                            '-r1', self.samples[sid]['r1'],
                                            '-r2', self.samples[sid]['r2'],
                                            '-avgq', '20'],
                          dependencies=['gunzip',])
        
@click.command()
@click.argument('run_dp')
def pre_analysis(run_dp):
    """ Pre-Analysis Management

    Sets up Pyflow WorkflowRunner with randomized log output string to avoid 
    possible output collision if ran multiple times.

    Arguments:
    run_dp -- String path to run directory to pre-analyze
    """
    log_output_dp = os.path.join(run_dp, 'bioinfo', "preanalysis")

    preanalyzer = PreAnalyzer(run_dp=run_dp)
    preanalyzer.run(mode='local', dataDirRoot=log_output_dp)


if __name__ == "__main__":
    pre_analysis()

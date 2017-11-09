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
    def __init__(self, analysis_dp, num_cpu=1):
        self.analysis_dp = analysis_dp
        self.read_dp = os.path.join(analysis_dp, 'reads')
        if not os.path.exists(self.read_dp):
            sys.exit('[ERROR]: No analysis reads, preanalysis must have failed!')
        self.num_cpu=num_cpu

    def generate_stability_file(self):
        fastqs = {}
        for fastq in os.listdir(self.read_dp):
            if '.fastq' not in fastq: continue
            sid = '_'.join([x for x in fastq.replace('-','_').split('_')][:-3])
            fastq = os.path.join(self.read_dp, fastq)
            
            if sid not in fastqs.keys():
                if '_R1_' in fastq: fastqs[sid] = {'r1': fastq}
                elif '_R2_' in fastq: fastqs[sid] = {'r2' :fastq}
                else: sys.exit('Fastq file {} neither for or rev!'.format(fastq))
            else:
                if '_R1_' in fastq: fastqs[sid]['r1'] = fastq
                elif '_R2_' in fastq: fastqs[sid]['r2'] = fastq
                else: sys.exit('Fastq file {} neither for or rev!'.format(fastq))
        output = open(os.path.join(self.read_dp, 'stability.files'), 'w+')
        for i, sid in enumerate(fastqs.keys()):
            if i != 0: output.write('\n')
            output.write('{}\t{}\t{}'.format(sid, fastqs[sid]['r1'], fastqs[sid]['r2']))
        output.close()
            

    def workflow(self):
        """ method invoked on class instance run call """
        self.generate_stability_file()
        
        # make contigs
        cmd = ['''mothur "#set.dir(input={});
                           make.contigs(file='stability.files',
                                        processors={});"'''.format(self.read_dp, self.num_cpu)]
        call(cmd, shell=True)

        cmd = ['''mothur "#set.dir(input={});
                           screen.seqs(fasta='stability.trim.contigs.fasta',
                                       group='stability.contigs.groups',
                                       maxambig=0,
                                       maxlength=275,
                                       minlength=240,
                                       maxhomop=8);"'''.format(self.read_dp)]
        call(cmd, shell=True)

        #m.unique.seqs(fasta='current')
        cmd = ['''mothur "#set.dir(input={});
                           unique.seqs(fasta='stability.trim.contigs.good.fasta');"'''.format(self.read_dp)]
        call(cmd, shell=True)

        #m.count.seqs(name='current', group='current')
        #m.align.seqs(fasta='current', reference='silva.seed_v123.pcr.align')
        #m.screen.seqs(fasta='current', count='current', start=1968, end=11550)
        cmd = ['''mothur "#set.dir(input={});
                  count.seqs(name='stability.trim.contigs.good.names',
                             group='stability.contigs.good.groups');
                  align.seqs(fasta='stability.trim.contigs.good.unique.fasta',
                             reference='silva.seed_v123.pcr.align', flip=T);
                  screen.seqs(fasta=current,
                              count=current,
                              optimize=start-end,
                              criteria=90);"'''.format(self.read_dp)]
        call(cmd, shell=True)

        #m.summary.seqs(fasta='current', count='current')
        cmd = ['''mothur "#set.dir(input={});
                           summary.seqs(fasta='stability.trim.contigs.good.unique.align');"'''.format(self.read_dp)]
        call(cmd, shell=True)

        cmd = ['''mothur "#set.dir(input={});
                           filter.seqs(fasta='stability.trim.contigs.good.unique.good.align',
                                       vertical=T, trump=.);"'''.format(self.read_dp)]
        call(cmd, shell=True)

        cmd = ['''mothur "#set.dir(input={});
                           unique.seqs(fasta='stability.trim.contigs.good.unique.good.filter.fasta',
                                       count='stability.trim.contigs.good.good.count_table');"'''.format(self.read_dp)]
        call(cmd, shell=True)

        cmd = ['''mothur "#set.dir(input={});
                           pre.cluster(fasta='stability.trim.contigs.good.unique.good.filter.unique.fasta',
                                       count='stability.trim.contigs.good.unique.good.filter.count_table',
                                       diffs=2);"'''.format(self.read_dp)]
        call(cmd, shell=True)

        # search for chimeras and remove them
        cmd = ['''mothur "#set.dir(input={});
                           chimera.uchime(fasta='stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta',
                                          count='stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table',
                                          dereplicate=t);
                           remove.seqs(fasta=current, accnos=current);
                           summary.seqs(count=current);"'''.format(self.read_dp)]
        call(cmd, shell=True)

        # classify sequences with a bayesian classifier using the silva database
        cmd = ['''mothur "#set.dir(input={});
                           classify.seqs(fasta='stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta',
                              count='stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table',
                                         reference='/nfs/vdenef-lab/Shared/Ruben/databases_taxass/silva.nr_v123.align',
                                         taxonomy=/nfs/vdenef-lab/Shared/Ruben/databases_taxass/silva.nr_v123.tax,
                                         cutoff=60);
                           remove.lineage(fasta=current,
                                          count=current,
                                          taxonomy=current,
                                          taxon=Chloroplast-Mitochondria-unknown-Eukaryota);
                           summary.seqs(count=current);
                           cluster.split(fasta=current,
                                         count=current,
                                         taxonomy=current,
                                         splitmethod=classify,
                                         taxlevel=4,
                                         cutoff=0.15);
                           make.shared(list=current,
                                       count=current,
                                       label=0.03);
                           classify.otu(list=current,
                                        count=current,
                                        taxonomy=current,
                                        label=0.03);"'''.format(self.read_dp)]
        call(cmd, shell=True)

        
        


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

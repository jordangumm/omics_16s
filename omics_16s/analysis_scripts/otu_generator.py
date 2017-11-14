import sys, os
import click

from subprocess import call


def generate_stability_file(read_dp):
    fastqs = {}
    for fastq in os.listdir(read_dp):
        if '.fastq' not in fastq: continue
        sid = '_'.join([x for x in fastq.replace('-','_').split('_')][:-3])
        fastq = os.path.join(read_dp, fastq)

        if sid not in fastqs.keys():
            if '_R1_' in fastq: fastqs[sid] = {'r1': fastq}
            elif '_R2_' in fastq: fastqs[sid] = {'r2' :fastq}
            else: sys.exit('Fastq file {} neither for or rev!'.format(fastq))
        else:
            if '_R1_' in fastq: fastqs[sid]['r1'] = fastq
            elif '_R2_' in fastq: fastqs[sid]['r2'] = fastq
            else: sys.exit('Fastq file {} neither for or rev!'.format(fastq))
    output = open(os.path.join(read_dp, 'stability.files'), 'w+')
    for i, sid in enumerate(fastqs.keys()):
        if i != 0: output.write('\n')
        output.write('{}\t{}\t{}'.format(sid, fastqs[sid]['r1'], fastqs[sid]['r2']))
    output.close()
    return None


@click.command()
@click.argument('read_dp')
@click.argument('dependencies_dp')
@click.option('--num_cpu', '-n', default=2)
def run(read_dp, dependencies_dp, num_cpu):
    """
    TODO: use built dependency silva files
    TODO: save outputs to independent directory for easier use by downstream processes
    """
    if not os.path.exists(read_dp):
        sys.exit('read_dp {} DOES NOT EXIST'.format(read_dp))
    if not os.path.exists(dependencies_dp):
        sys.exit('dependencies_dp {} DOES NOT EXIST'.format(dependencies_dp))
    generate_stability_file(read_dp)

    if not os.path.exists(os.path.join(dependencies_dp, 'silva', 'silva.seed_v128.pcr.align')):
        cmd = ['''mothur "#pcr.seqs(fasta={}, start=11894, end=25319, keepdots=F);"'''.format(
                             os.path.join(dependencies_dp, 'silva', 'silva.seed_v128.align'))]
        call(cmd, shell=True) 

    # make contigs
    cmd = ['''mothur "#set.dir(input={});
                       make.contigs(file='stability.files',
                                    processors={});
                       screen.seqs(fasta=current,
                                   group=current,
                                   maxambig=0,
                                   maxlength=275,
                                   minlength=240,
                                   maxhomop=8);
                       summary.seqs(count=current);
                       unique.seqs(fasta=current);
                       count.seqs(name=current,
                                  group=current);
                       align.seqs(fasta=current,
                                  reference={});
                       summary.seqs(count=current);
                       screen.seqs(fasta=current,
                                   count=current,
                                   start=1968,
                                   end=11550);
                       summary.seqs(count=current);
                       filter.seqs(fasta=current,
                                   vertical=T,
                                   trump=.);
                       unique.seqs(fasta=current,
                                   count=current);
                       pre.cluster(fasta=current,
                                   count=current,
                                   diffs=2);
                       chimera.uchime(fasta=current,
                                      count=current,
                                      dereplicate=t);
                       remove.seqs(fasta=current, accnos=current);
                       summary.seqs(count=current);
                       classify.seqs(fasta=current,
                                     count=current,
                                     reference={},
                                     taxonomy={},
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
                                    label=0.03);"'''.format(read_dp, num_cpu,
              os.path.join(dependencies_dp, 'silva', 'silva.seed_v128.pcr.align'),
              os.path.join(dependencies_dp, 'silva', 'silva.nr_v128.align'),
              os.path.join(dependencies_dp, 'silva', 'silva.nr_v128.tax'))]
    call(cmd, shell=True)
    return None


if __name__ == "__main__":
    run()

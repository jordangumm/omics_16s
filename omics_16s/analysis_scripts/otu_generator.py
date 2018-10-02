import sys, os
import screed
import click

from subprocess import call


def generate_stability_file(read_dp, output_dp):
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
    output = open(os.path.join(output_dp, 'stability.files'), 'w+')
    for i, sid in enumerate(fastqs.keys()):
        if i != 0: output.write('\n')
        output.write('{}\t{}\t{}'.format(sid, fastqs[sid]['r1'], fastqs[sid]['r2']))
    output.close()
    return None


@click.command()
@click.argument('reads_dp')
@click.argument('mothur_dp')
@click.argument('dependencies_dp')
@click.option('--num_cpu', '-n', default=2)
def run(reads_dp, mothur_dp, dependencies_dp, num_cpu):
    """
    TODO: save outputs to independent directory for easier use by downstream processes
    """
    if not os.path.exists(reads_dp):
        sys.exit('read_dp {} DOES NOT EXIST'.format(reads_dp))
    if not os.path.exists(mothur_dp):
        os.makedirs(mothur_dp)
    if not os.path.exists(dependencies_dp):
        sys.exit('dependencies_dp {} DOES NOT EXIST'.format(dependencies_dp))
    generate_stability_file(reads_dp, mothur_dp)

    if not os.path.exists(os.path.join(dependencies_dp, 'silva', 'silva.seed_v132.pcr.align')):
        cmd = ['''mothur "#pcr.seqs(fasta={}, start=11894, end=25319, keepdots=F, processors={});"'''.format(
                             os.path.join(dependencies_dp, 'silva', 'silva.seed_v132.align'), num_cpu)]
        call(cmd, shell=True) 
    if not os.path.exists(os.path.join(dependencies_dp, 'silva', 'silva.seed_v132.pcr.align')):
        sys.exit('silva.seed_v132.pcr.align NOT CREATED.')

    miniconda_bin_dp = os.path.join(dependencies_dp, 'mothur', 'bin')

    # make contigs
    cmd = ['''mothur "#set.dir(input={}, tempdefault={});
                       make.contigs(file='stability.files',
                                    processors={});
                       screen.seqs(fasta=current,
                                   group=current,
                                   maxambig=0,
                                   maxlength=275,
                                   minlength=240);
                       summary.seqs(count=current);
                       unique.seqs(fasta=current);
                       count.seqs(name=current,
                                  group=current);
                       align.seqs(fasta=current,
                                  reference={});
                       summary.seqs(fasta=current,
                                    count=current);
                       screen.seqs(fasta=current,
                                   count=current,
                                   summary=current,
                                   start=1968,
                                   end=11550,
                                   maxhomop=8);
                       summary.seqs(count=current,
                                    count=current);
                       filter.seqs(fasta=current,
                                   vertical=T,
                                   trump=.);
                       unique.seqs(fasta=current,
                                   count=current);
                       pre.cluster(fasta=current,
                                   count=current,
                                   diffs=2);
                       chimera.vsearch(fasta=current,
                                      count=current,
                                      dereplicate=t);
                       remove.seqs(fasta=current, accnos=current);
                       summary.seqs(fasta=current,
                                    count=current);
                       classify.seqs(fasta=current,
                                     count=current,
                                     reference={},
                                     taxonomy={},
                                     cutoff=80);
                       remove.lineage(fasta=current,
                                      count=current,
                                      taxonomy=current,
                                      taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);
                       summary.seqs(fasta=current,
                                    count=current);
                       cluster.split(fasta=current,
                                     count=current,
                                     taxonomy=current,
                                     splitmethod=classify,
                                     taxlevel=4,
                                     cutoff=0.03);
                       make.shared(list=current,
                                   count=current,
                                   label=0.03);
                       classify.otu(list=current,
                                    count=current,
                                    taxonomy=current,
                                    label=0.03);
                       tree.shared(shared=current,
                                   calc=jest-thetayc-braycurtis);
                       get.oturep(column=current,
                                  name=current,
                                  fasta=current,
                                  list=current);"'''.format(mothur_dp, miniconda_bin_dp, num_cpu,
              os.path.join(dependencies_dp, 'silva', 'silva.seed_v132.pcr.align'),
              os.path.join(dependencies_dp, 'silva', 'silva.nr_v132.align'),
              os.path.join(dependencies_dp, 'silva', 'silva.nr_v132.tax'))]
    call(cmd, shell=True)
    #db = screed.read_fasta_sequences(os.path.join(mothur_dp, 'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.opti_mcc.unique_list.0.03.rep.fasta'))
    db = screed.read_fasta_sequences(os.path.join(mothur_dp, 'stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.rep.fasta'))
    output = open(os.path.join(mothur_dp, '../', 'otus.fasta'), 'w+')
    for otu in db:
        output.write('>{}\n'.format(otu.split('\t')[1].split('|')[0]))
        output.write('{}\n'.format(db[otu].sequence))
    output.close()
    return None


if __name__ == "__main__":
    run()

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
    def __init__(self, analysis_dp, taxass_dp, silva_dp, num_cpu):
        self.output_dp = os.path.join(analysis_dp, 'taxass')
        if not os.path.exists(self.output_dp):
            os.makedirs(self.output_dp)

        self.fasta_fp = os.path.join(analysis_dp, 'mothur',
             'stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta')
        self.count_fp = os.path.join(analysis_dp, 'mothur',
             'stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table')
        self.names_fp = os.path.join(analysis_dp, 'mothur',
             'stability.trim.contigs.good.names')
        if not os.path.exists(self.fasta_fp) or not os.path.exists(self.count_fp) or not os.path.exists(self.names_fp):
            sys.exit('mothur files not available for taxass analysis')

        self.taxass_dp = taxass_dp # taxass scripts
        self.specific_fasta = os.path.join(taxass_dp, 'FreshTrain18Aug2016', 'FreshTrain18Aug2016.fasta')
        self.specific_taxa = os.path.join(taxass_dp, 'FreshTrain18Aug2016', 'FreshTrain18Aug2016.taxonomy')
        if not os.path.exists(self.specific_fasta) or not os.path.exists(self.specific_taxa):
            sys.exit('[ERROR]: Taxass file not found in {}'.format(taxass_dp))

        self.calc_full_length_fp = os.path.join(taxass_dp, 'tax-scripts', 'calc_full_length_pident.R')
        self.filter_seqids_fp = os.path.join(taxass_dp, 'tax-scripts', 'filter_seqIDs_by_pident.R')
        self.plot_blast_hits_fp = os.path.join(taxass_dp, 'tax-scripts', 'plot_blast_hit_stats.R')
        self.find_seqids_fp = os.path.join(taxass_dp, 'tax-scripts', 'find_seqIDs_blast_removed.py')
        self.create_fastas_fp = os.path.join(taxass_dp, 'tax-scripts', 'create_fastas_given_seqIDs.py')

        self.general_fasta = os.path.join(silva_dp, 'silva.nr_v128.align')
        self.general_taxa = os.path.join(silva_dp, 'silva.nr_v128.tax')
        self.num_cpu=num_cpu

        self.taxass_fasta = os.path.join(self.output_dp, 'sequence.fasta')


    def workflow(self):
        """ method invoked on class instance run call """
        cmd = "sed -e '/>/!s/-//g' < $(echo {}) > {}".format(self.fasta_fp, self.taxass_fasta)
        self.addTask('make_taxass_fasta', command=cmd)

        custom_db = os.path.join(self.output_dp, 'custom_database', 'FWonly_18Aug2016custom.db')
        cmd = "makeblastdb -dbtype nucl -in $(echo {}) -input_type fasta -parse_seqids -out {}".format(self.specific_fasta, custom_db)
        self.addTask('make_customdb', command=cmd, dependencies=['make_taxass_fasta',])
        
        custom_blast = os.path.join(self.output_dp, 'custom.blast')
        cmd = "blastn -query {} -task megablast -db {} -out {} -outfmt 11 -max_target_seqs 5 -num_threads {}".format(self.taxass_fasta, custom_db, custom_blast, self.num_cpu)
        self.addTask('blastn', command=cmd, dependencies=['make_customdb',])

        blast_table = os.path.join(self.output_dp, 'otus.custom.blast.table')
        cmd = 'blast_formatter -archive {} -outfmt "6 qseqid pident length qlen qstart qend" -out {}'.format(custom_blast, blast_table)
        self.addTask('blast_formatter', command=cmd, dependencies=['blastn',])

        blast_table_modified = os.path.join(self.output_dp, 'otus.custom.blast.table.modified')
        cmd = 'Rscript $(echo "{}") {} {}'.format(self.calc_full_length_fp, blast_table, blast_table_modified)
        self.addTask('correct_blast', command=cmd, dependencies=['blast_formatter',])

        ids_above_fp = os.path.join(self.output_dp, 'ids.above.97')
        ids_below_fp = os.path.join(self.output_dp, 'ids.below.97')
        cmd = 'Rscript $(echo "{}") {} {} 97 TRUE'.format(self.filter_seqids_fp, blast_table_modified, ids_above_fp)
        self.addTask('filter_above', command=cmd, dependencies=['correct_blast',])
        cmd = 'Rscript $(echo "{}") {} {} 97 FALSE'.format(self.filter_seqids_fp, blast_table_modified, ids_below_fp)
        self.addTask('filter_below', command=cmd, dependencies=['correct_blast',])

        cmd = 'Rscript $(echo "{}") {} 97 {}'.format(self.plot_blast_hits_fp, blast_table_modified, self.output_dp)
        self.addTask('make_plots', command=cmd, dependencies=['filter_above', 'filter_below',])

        ids_missing_fp = os.path.join(self.output_dp, 'ids.missing')
        cmd = 'python $(echo "{}") {} {} {}'.format(self.find_seqids_fp, self.taxass_fasta, blast_table_modified, ids_missing_fp)
        self.addTask('calc_ids_missing', command=cmd, dependencies=['make_plots',])

        ids_below_all_fp = os.path.join(self.output_dp, 'ids.below.97.all')
        cmd = 'cat {} {} > {}'.format(ids_below_fp, ids_missing_fp, ids_below_all_fp)
        self.addTask('cat_ids_below', command=cmd, dependencies=['calc_ids_missing',])

        above_final_fp = os.path.join(self.output_dp, 'otus.above.97.fasta')
        below_final_fp = os.path.join(self.output_dp, 'otus.below.97.fasta')
        cmd = 'python $(echo "{}") {} {} {}'.format(self.create_fastas_fp, ids_above_fp, self.taxass_fasta, above_final_fp)
        self.addTask('create_above_final', command=cmd, dependencies=['filter_above',])

        cmd = 'python $(echo "{}") {} {} {}'.format(self.create_fastas_fp, ids_below_all_fp, self.taxass_fasta, below_final_fp)
        self.addTask('create_below_final', command=cmd, dependencies=['cat_ids_below',])

        cmd= 'mothur "#classify.seqs(seed=777, fasta={}, template={},  taxonomy={}, method=wang, probs=T, processors={}, cutoff=80)"'.format(above_final_fp, self.specific_fasta, self.specific_taxa, self.num_cpu)
        self.addTask('assign_above_taxonomy', command=cmd, dependencies=['create_above_final',])

        cmd= 'mothur "#classify.seqs(seed=777, fasta={}, template={},  taxonomy={}, method=wang, probs=T, processors={}, cutoff=80)"'.format(below_final_fp, self.general_fasta, self.general_taxa, self.num_cpu)
        self.addTask('assign_below_taxonomy', command=cmd, dependencies=['create_below_final',])

        # created in assign_above_taxonomy
        above_taxonomy_fp = os.path.join(self.output_dp, 'otus.above.97.FreshTrain18Aug2016.wang.taxonomy')
        below_taxonomy_fp = os.path.join(self.output_dp, 'otus.below.97.nr_v128.wang.taxonomy')
        combined_taxonomy_fp = os.path.join(self.output_dp, 'final.FWDB.Silva.taxonomy')
        cmd = 'cat {} {} > {}'.format(above_taxonomy_fp, below_taxonomy_fp, combined_taxonomy_fp)
        self.addTask('combine_taxonomy', command=cmd, dependencies=['assign_above_taxonomy', 'assign_below_taxonomy'])



@click.command()
@click.argument('analysis_dp')
@click.argument('taxass_dp')
@click.argument('silva_dp')
@click.option('--num_cpu', '-n', type=click.INT, default=2)
def analysis(analysis_dp, taxass_dp, silva_dp, num_cpu):
    """ Runs mothur fasta output against TaxAss workflow

    \b
    Arguments:
        analysis_dp: directory path to active analysis directory
        fasta_fp: file path to mothur aligned fasta file
        taxass_dp: directory path to taxass software install
        silva_dp: directory path to mothur released silva files
    """
    log_output_dp = os.path.join(analysis_dp, 'logs', 'taxass')

    taxasser = TaxAss(analysis_dp=analysis_dp,
                           taxass_dp=taxass_dp,
                           silva_dp=silva_dp,
                           num_cpu=num_cpu)
    taxasser.run(mode='local', dataDirRoot=log_output_dp)

      
    combined_taxonomy_fp = os.path.join(taxasser.output_dp, 'final.FWDB.Silva.taxonomy')
    cmd = ['''mothur "#set.dir(output={});
                       remove.lineage(fasta={},
                                      count={},
                                      taxonomy={},
                                      taxon=Mitochondria-unknown-Archaea-Eukaryota-Chloroplast);
                        cluster.split(fasta=current,
                                      count=current,
                                      taxonomy=current,
                                      splitmethod=classify,
                                      taxlevel=4,
                                      cutoff=0.03,
                                      method=average);
                        summary.seqs(count=current);
                        make.shared(list=current,
                                    count=current,
                                    label=0.03);
                        classify.otu(list=current,
                                     count=current,
                                     taxonomy={},
                                     label=0.03);
                        get.lineage(taxonomy={},
                                    taxon=Bacteria,
                                    name={});
                        list.seqs(taxonomy=current);
                        get.seqs(accnos=current,
                                 fasta={},
                                 name=current);
                        deunique.seqs(fasta=current,
                                      name=current);
                        get.oturep(list=current,
                                   method=abundance,
                                   fasta={},
                                   name={},
                                   label=0.03);
                        tree.shared(shared=current,
                             calc=jest-thetayc-braycurtis);"'''.format(taxasser.output_dp,
                                                         taxasser.fasta_fp,
                                                         taxasser.count_fp,
                                                         combined_taxonomy_fp,
                                                         combined_taxonomy_fp,
                                                         combined_taxonomy_fp,
                                                         taxasser.names_fp,
                                                         taxasser.fasta_fp.replace('.precluster.pick.pick',''),
                                                         taxasser.fasta_fp,
                                                         taxasser.names_fp)]
    call(cmd, shell=True)
    return None



if __name__ == "__main__":
    analysis()

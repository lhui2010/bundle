"""
Various assembler wrapper
"""
import os

import os.path as op

from iga.apps.base import conda_act, Config, mkdir, get_prefix, sh, bsub, emain, abspath_list, waitjob, logger


def bam2fastq(subreads=None):
    """
    Input zeins.bam zeins.bam.pbi
    Output zeins.fasta.gz
    :param subreads:
    :return:
    """
    if type(subreads) == list:
        abspath_list(subreads)
        subreads = " ".join(subreads)
    else:
        subreads = op.abspath(subreads)

    if not op.exists(subreads + '.pbi'):
        logger.error("{}.pbi is needed! Exiting..".format(subreads))
        return 1

    prefix = get_prefix(subreads)
    cmd = 'bam2fasta -o {0} {1}'.format(prefix, subreads)
    job = bsub(cmd, name='bam2fasta')
    waitjob(job)
    return '{0}.fasta.gz'.format(prefix)


# 0 cfg_file
falcon_sh = r"""
fc_run {} >>falcon_run.out 2>>falcon_run.err
"""


def falcon(subreads=None, genome_size=None, prefix='', etc=''):
    r"""
    :param subreads: pacbio DNA subreads (bam accepted)
    :param genome_size: (genomesize in bp)
    :param prefix: (species name)
    :param etc: (other fields need to be updated in falcon cfg)
    :return:
    """
    logger.debug(subreads)

    if '.bam' in subreads:
        subreads = bam2fastq(subreads)

    logger.debug(subreads)

    if type(subreads) == list:
        abspath_list(subreads)
        subreads = " ".join(subreads)
    else:
        subreads = op.abspath(subreads)

    if prefix == '':
        prefix = get_prefix(subreads)

    workdir = 'workdir_falcon_{}'.format(prefix)
    mkdir(workdir)
    os.chdir(workdir)

    fofn_file = 'bam.fofn'
    cfg_file = "falcon.cfg"
    fofn_file = op.abspath(fofn_file)
    sh("ls {0} > {1}".format(subreads, fofn_file))
    cfg = Config('falcon')
    cfg.update('[General]input_fofn={0}'.format(fofn_file))
    cfg.update('[General]genome_size={0}'.format(genome_size))
    if etc != '':
        cfg.update(etc)
    cfg.write_to_file(cfg_file)

    cmd = conda_act.format('falcon') + falcon_sh.format(cfg_file)
    bsub(cmd, name='falcon' + prefix)


# 0 reference contig (abs path)
# 1 subreads.fasta (abs path)
purge_dups_sh = r"""
pri_asm={0}
subreads={1}


minimap2 -t 40 -xmap-pb $pri_asm $subreads | gzip -c - > align.paf.gz

pbcstat *.paf.gz #(produces PB.base.cov and PB.stat files)
calcuts PB.stat > cutoffs 2>calcults.log
#        Notice If you have a large genome, please set minimap2 -I option to ensure the genome can be 
# indexed once, otherwise read depth can be wrong.

#        Step 1. Split an assembly and do a self-self alignment. Commands are following:
split_fa $pri_asm > $pri_asm.split
minimap2 -t 40 -xasm5 -DP $pri_asm.split $pri_asm.split | gzip -c - > $pri_asm.split.self.paf.gz

#Step 2. Purge haplotigs and overlaps with the following command.
purge_dups -2 -T cutoffs -c PB.base.cov $pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log

#Step 3. Get purged primary and haplotig sequences from draft assembly.
get_seqs dups.bed $pri_asm
"""


def purge_dups(contig=None, subreads=None, prefix=''):
    r"""
    purge dups wrapper
    :param contig:  contig fasta
    :param subreads: subreads.fasta[.gz]
    :param prefix: usually name of contig.fasta
    :return:
    """
    contig = op.abspath(contig)
    subreads = op.abspath(subreads)
    if prefix == '':
        prefix = get_prefix(subreads)

    workdir = 'workdir_purge_dups_{}'.format(prefix)
    mkdir(workdir)
    os.chdir(workdir)

    cmd = purge_dups_sh.format(contig, subreads)
    bsub(cmd, name="purge_dups_" + prefix)
    return 0


if __name__ == '__main__':
    emain()

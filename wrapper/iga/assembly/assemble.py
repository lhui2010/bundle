"""
Various assembler wrapper
"""
import os

import os.path as op

from iga.apps.base import conda_act, Config, mkdir, get_prefix, sh, bsub, emain

# 0 subreads.fasta
# 1 prefix
# 2 fofn file
falcon_sh = r"""
fc_run ${cfg} >>falcon_run.out 2>>falcon_run.err
"""


def falcon(subreads=None, genome_size=None, prefix='', etc=''):
    r"""
    :param subreads: pacbio DNA subreads (bam accepted)
    :param genome_size: (genomesize in bp)
    :param prefix: (species name)
    :param etc: (other fields need to be updated in falcon cfg)
    :return:
    """
    if type(subreads) == list:
        subreads = " ".join(subreads)
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
    cfg.update('input_fofn={0}'.format(fofn_file))
    cfg.update('genome_size={0}'.format(genome_size))
    if etc != '':
        cfg.update(etc)
    cfg.write_to_file(cfg_file)

    cmd = conda_act.format('falcon') + falcon_sh.format(subreads, cfg_file)
    bsub(cmd, name='falcon' + prefix)


if __name__ == '__main__':
    emain()

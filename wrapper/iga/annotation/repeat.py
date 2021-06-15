"""
Repeat annotation utils
"""
import os

import os.path as op

from iga.apps.base import emain, conda_act, bsub

# def repeat_mask():


def goto_workdir(program, sample=''):
    """
    :param program:
    :param sample:
    :return:
    """
    folder_name = 'workdir_{}_{}'.format(program, sample)
    cmd = "mkdir -p {0}; cd {0}\n".format(folder_name)
    return cmd


def edta(genome=None, threads=20):
    """
    :param genome:
    :return:
    """
    abs_genome = op.abspath(genome)
    rel_genome = op.basename(genome)
    cmd = conda_act.format('EDTA')
    cmd += goto_workdir('edta', op.basename(rel_genome))
    cmd += 'EDTA.pl --anno 1  --genome {0} --threads {1} \
    >edta_anno.out 2>edta_anno.err'.format(abs_genome, threads)
    bsub(cmd, name='EDTA{}'.format(rel_genome), cpus=threads)


if __name__ == "__main__":
    emain()

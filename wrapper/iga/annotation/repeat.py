"""
Repeat annotation utils
"""
import logging
import os

import os.path as op

from iga.apps.base import emain, conda_act, bsub

# def repeat_mask():


# 0 ref.fa
# 1 repeat.gff
# 2 threads
repeat_masker_sh = r"""
BuildDatabase -name {0} -engine ncbi {0} 
RepeatModeler -engine ncbi -pa {1} -database {0}
mkdir -p {0}.custom_lib.out
RepeatMasker -lib ref-families.fa {0} -pa {1} -dir {0}.custom_lib.out"
"""


# species Viridiplantae
def repeatmasker(genome=None, species='', denovo='T', threads=30):
    """
    :param genome:
    :param species:
    :param denovo:
    :return:
    """

    cmd = goto_workdir('repeatmask', sample=genome)
    genome = "../{}".format(genome)
    if species != '':
        cmd += "\nmkdir -p species_lib.out && RepeatMasker {0} -species {1} -pa {2} -dir species_lib.out\n".format(
            genome, species, threads
        )
    elif denovo == 'T':
        cmd += repeat_masker_sh.format(genome, threads)
    else:
        logging.error("Either provide a species name or use denovo prediction mode")
        exit(1)
    bsub(cmd, name="repeat_masker_{}".format(genome), threads=threads)


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

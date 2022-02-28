"""
gemoma relevant utils
"""
import logging
from os.path import abspath

from iga.annotation.maker import maker_rename_gff
from iga.annotation.repeat import goto_workdir
from iga.apps.base import conda_act, bsub, waitjob, get_prefix, emain, abspath_list, sh

# 0 input_genome
# 1 ortho_genome
# 2 ortho_gff
# 3 threads
# 4 output
gemoma_pipe_sh = r"""
GeMoMa GeMoMaPipeline threads={3} outdir={4} GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t={0} i={4} a={2} g={1}


# GeMoMa GeMoMaPipeline threads=128 outdir=Mt.outdir GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=Duparquetia_orchidacea.fa.masked.fa i=Mt a=Medicago_truncatula.gff g=Medicago_truncatula.fa"
# GeMoMa GeMoMaPipeline threads=128 outdir=At.outdir GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=Duparquetia_orchidacea.fa.masked.fa i=At a=Arabidopsis_thaliana.gff g=Arabidopsis_thaliana.fa"
"""


def gemoma_run(input_genome=None, ortho_genome=None, ortho_gff=None, output='', threads=40):
    """
    The GeMoMaPipeline wrapper (assuming you have installed gemoma in conda env gemoma)
    Args:
        input_genome:
        ortho_genome:
        ortho_gff:
        output:
        threads:
    Returns:
    """
    if output == '':
        output = get_prefix(input_genome) + '.' + get_prefix(ortho_genome)
    cmd = conda_act.format('gemoma')
    cmd += gemoma_pipe_sh.format(input_genome, ortho_genome, ortho_gff, threads, output)
    jobid = bsub(cmd, name="gemoma_{}".format(output), cpus=threads)
    return jobid


#gemoma_gaf_sh = r"""
#GeMoMa GAF g={0}/final_annotation.gff g={1}/final_annotation.gff g={2}/final_annotation.gff
#"""


def gemoma_gaf(dirs=None, outdir=''):
    cmd = conda_act.format('gemoma')
    cmd += "GeMoMa GAF "
    for d in dirs:
        cmd += 'g={}/final_annotation.gff '.format(d)
    jobid = bsub(cmd, name='gemoma_gaf')
    waitjob(jobid)
    logging.info('Merged gff file is: filtered_predictions.gff')


def gemoma_pipe(input_genome=None, homo_genome='', homo_gff='', threads=40, outdir=''):
    """
    This is a wrapper that combines GeMoMaPipeline, GAF, and gff rename
    homo_genome: seperated with ','
    Returns:
    """
    # cmd = goto_workdir(input_genome+'.gemoma')
    input_prefix = get_prefix(input_genome)

    homo_genomes = homo_genome.split(',')
    homo_gffs = homo_gff.split(',')

    # sh(cmd)
    # for g in homo_genomes:
    #     sh('ln -s ../{}'.format(g))

    jobids = []
    #for i in range(0, len(homo_genomes)):
    #    jobid = gemoma_run(input_genome=input_genome, ortho_genome=homo_genomes[i], ortho_gff=homo_gffs[i], output='', threads=threads)
    #    jobids.append(jobid)
    #waitjob(jobids)

    outdirs = []
    for i in homo_genomes:
        outdirs.append("{}.{}".format(get_prefix(input_genome), get_prefix(i)))
    gemoma_gaf(dirs=outdirs)

    raw_gff = input_prefix + ".gff"
    format_gff = input_prefix + ".format.gff"
    sh('mv filtered_predictions.gff {}'.format(raw_gff))
    sh('mv protocol_GAF.txt {}.protocol_GAF.txt'.format(input_prefix))

    (genus, species) = input_genome.split('.')[0].split('_')
    gff_prefix = genus[:2] + species[:3]
    gff_prefix = gff_prefix.upper()
    maker_rename_gff(gff=raw_gff, prefix=gff_prefix)
    logging.info("The resulting gff is {}".format(format_gff))
    pass


if __name__ == '__main__':
    emain()

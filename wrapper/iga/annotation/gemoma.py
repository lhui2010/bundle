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
GeMoMa GeMoMaPipeline threads={3} outdir={4} GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t={0} i={4} a={2} g={1}"
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
        output = get_prefix(ortho_genome)
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
    homo_genomes = homo_genome.split(',')
    homo_gff = homo_gff.split(',')
    homo_genomes = abspath_list(homo_genomes)
    homo_gff = abspath_list(homo_gff)
    cmd = goto_workdir(input_genome+'.gemoma')
    sh(cmd)
    jobids = []
    for i in range(0, len(homo_genomes)):
        jobid = gemoma_run(input_genome=input_genome, ortho_genome=homo_genomes[i], ortho_gff=homo_gff[i], output='', threads=threads)
        jobids.append(jobid)
    waitjob(jobids)

    outdirs = []
    for i in homo_genomes:
        outdirs.append(get_prefix(i))
    gemoma_gaf(dirs=outdirs)

    (genus, species) = input_genome.split('.')[0].split('_')
    gff_prefix = genus[:2] + species[:3]
    gff_prefix = gff_prefix.upper()
    maker_rename_gff(gff='filtered_predictions.gff', prefix=gff_prefix)
    sh("cp filtered_predictions.gff ../{0}.gemoma.gff")
    pass


if __name__ == '__main__':
    emain()

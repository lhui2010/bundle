"""
gemoma relevant utils
"""
import logging

from iga.annotation.maker import maker_rename_gff
from iga.apps.base import conda_act, bsub, waitjob, get_prefix

# 0 input_genome
# 1 ortho_genome
# 2 ortho_gff
# 3 threads
# 4 output
gemoma_pipe_sh = r"""
GeMoMa GeMoMaPipeline threads={3} outdir={4} GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t={0} i={4} a={2} g={1}"
GeMoMa GeMoMaPipeline threads=128 outdir=Mt.outdir GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=Duparquetia_orchidacea.fa.masked.fa i=Mt a=Medicago_truncatula.gff g=Medicago_truncatula.fa"
GeMoMa GeMoMaPipeline threads=128 outdir=At.outdir GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=Duparquetia_orchidacea.fa.masked.fa i=At a=Arabidopsis_thaliana.gff g=Arabidopsis_thaliana.fa"
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
    pass


gemoma_gaf_sh = r"""
GeMoMa GAF g=gemoma_outdir/final_annotation.gff g=Mt.outdir/final_annotation.gff g=At.outdir/final_annotation.gff
"""


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
    Returns:
    """
    gemoma_run(1)
    gemoma_run(2)
    gemoma_gaf()
    maker_rename_gff()
    pass

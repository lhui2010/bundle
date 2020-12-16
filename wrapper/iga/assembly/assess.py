#!/usr/bin/env python

"""
running busco assessement on genomes
"""

import argparse
import textwrap
import subprocess
import os
from iga.apps.base import emain, conda_act, bsub, sh, waitjob

# class CTL():

# Give arguments, will return the specific CTLs
# Support Direct Print
# Support tag value change


# 0 input fasta (relative dir)
# 2 threads
lai_sh = """
PREFIX=`basename {0}`
mkdir -p workdir_LAI_$PREFIX && cd  workdir_LAI_$PREFIX
ln -s {0}
LTR_FINDER_parallel -seq $PREFIX -threads {1} -harvest_out -size 1000000 -time 300
cat $PREFIX.harvest.scn $PREFIX.finder.combine.scn > $PREFIX.rawLTR.scn
LTR_retriever -genome $PREFIX -inharvest $PREFIX.rawLTR.scn -threads {1} 
LAI -genome $PREFIX -intact $PREFIX.pass.list -all $PREFIX.out {1}"
echo -n "LAI for {0} is: "
sed -n '2p' $PREFIX.out.LAI |awk '{print $7}'
"""

def lai(genome=None, threads=30):
    r"""
    Calculate lai for specific
    :param genome:
    :param threads:
    :return:
    """
    cmd = conda_act.format("EDTA")
    cmd += lai_sh.format(genome, threads)
    job = bsub(cmd)
    waitjob(job)
    return 0

# 0 threads
# 1 mode genome/pep
# 2 input
# 3 output
# 4 lineage

busco_sh = """
busco -f -c {0} -m {1} -i {2} -o {3} -l {4}
"""


def busco(genome_fasta=None, mode='genome', lineage='embryophyta_odb10', threads=64, output='', usegrid='T'):
    """
    :param genome_fasta:
    :param mode: genome | prot
    :param lineage:
    :param threads:
    :param output:
    :return:
    """
    busco_export = r"""
    #export AUGUSTUS_CONFIG_PATH=/tmp/lh_config
    export BUSCO_CONFIG_FILE=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/busco/myconfig.ini
    """
    # deploy_augustus = r"""
    # touch /tmp/lh_config && rm -fr /tmp/lh_config && cp -fr  /ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/config /tmp/lh_config
    # """
    if (output == ''):
        output = os.path.basename(genome_fasta) + ".busco.embryophyta.v4.1.2"
    #cmd = conda_act.format('busco') + deploy_augustus + busco_sh.format(threads, mode, genome_fasta, output, lineage)
    cmd = conda_act.format('busco') + busco_export + busco_sh.format(threads, mode, genome_fasta, output, lineage)
    #    cmd = conda_act
    if(usegrid == 'T'):
        bsub(cmd, cpus=5)
    else:
        sh(cmd)


#
# def main():
#     prog_name = "busco_wrapper"
#     usage = "run busco on selected GENOME"
#
#     parser = argparse.ArgumentParser(
#         prog=prog_name,
#         formatter_class=argparse.RawDescriptionHelpFormatter,
#         description=textwrap.dedent(usage),
#         epilog="")
#     parser.add_argument("GENOME", help="Genome to be evalutated in fasta format")
#     parser.add_argument("-t", "--threads", default=64, type=int, help="flanking distance default (1000)")
#     args = parser.parse_args()
#
#     busco(args.GENOME)
# #    flanking_distance = args.flanking

if __name__ == "__main__":
    emain()

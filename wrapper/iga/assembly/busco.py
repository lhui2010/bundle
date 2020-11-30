#!/usr/bin/env python

"""
running busco assessement on genomes
"""

import argparse
import textwrap
import subprocess
import os

# class CTL():

# Give arguments, will return the specific CTLs
# Support Direct Print
# Support tag value change


# 0 threads
# 1 mode genome/pep
# 2 input
# 3 output
# 4 lineage
from iga.apps.base import emain

busco_sh = """
busco -f -c {0} -m {1} -i {2} -o {3} -l {4}
"""


def busco(genome_fasta=None, mode='genome', lineage='embryophyta_odb10', threads=64, output=''):
    # assembly, subreads
    conda_act = r"""
    export AUGUSTUS_CONFIG_PATH=/tmp/lh_config
    export BUSCO_CONFIG_FILE=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/busco/myconfig.ini
    """
    deploy_augustus = r"""
    touch /tmp/lh_config && rm -fr /tmp/lh_config && cp -fr  /ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/config /tmp/lh_config
    """
    if (output == ''):
        output = os.path.basename(genome_fasta) + ".busco.embryophyta.v4.1.2"
    cmd = conda_act.format('busco') + deploy_augustus + busco_sh.format(threads, mode, genome_fasta, output, lineage)
    #    cmd = conda_act
    subprocess.run(cmd, shell=True)


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

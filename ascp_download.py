#!/usr/bin/env python
import argparse
import os
import logging
import sys




#Args
parser = argparse.ArgumentParser(description='SRR downloader with ascp')
parser.add_argument('SRR_list', type=str, nargs = 1,
                    help='file contain SRR ID in each line')
args = parser.parse_args()

outdir = "~/ncbi/public/sra"
outdir = "./"

def srr2ascp(srr_id):
    """Transform srr ID to download command """
    "ascp -T -l640M -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:1GB /tmp/"
    "ascp -T -l640M -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR212/SRR2120187/SRR2120187.sra  /tmp/"
    prefix = "ascp -T -l640M -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/"
    middle = "/".join([srr_id[:6], srr_id, srr_id + ".sra"])
    suffix = ' '+ outdir

    return prefix + middle + suffix
    
##Logging


#Prepare fasta
if __name__ == '__main__':
#https://stackoverflow.com/questions/14058453/making-python-loggers-output-all-messages-to-stdout-in-addition-to-log-file
#Format logging
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)

    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.WARNING)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    root.addHandler(handler)

    logging.basicConfig(format='%(asctime)s %(message)s')
    with open(args.SRR_list[0]) as fh:
        for line in fh:
            ascp_cmd = srr2ascp(line.rstrip())
            logging.warning("Downloading " + line.rstrip())
            os.system(ascp_cmd)
            logging.warning("Finished " + line.rstrip())

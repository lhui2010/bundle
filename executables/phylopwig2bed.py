import sys
import os.path

### from msa-pipeline

def phylopwig2bed(inwig,outbed):
    '''
    Convert phyloP base-by-base wig output to bed format.
    Checks if chrom matches input file name and if not
    replaces wig chrom with input file name chrom to 
    prevent phyloP stripping chrom name, e.g. after "." char.
    Bedops wig2bed alternative was not robust in the
    snakemake setting.
    '''
    full_chrom = os.path.basename(inwig)[:-11]
    with open(outbed,'a') as out:
        with open(inwig,'r') as w:
            # definition line
            header = next(w)
            # split on whitespace
            header = header.split()
            # strip leading "chrom="
            chrom = header[1][6:]
            if chrom != full_chrom:
                chrom = full_chrom
            # strip leading "start="
            # and convert 1-indexed to 0-indexed BED position
            bedstop = int(header[2][6:])
            bedstart = bedstop - 1
            for line in w:
                line = line.strip()
                if line.startswith("fixedStep"):
                    newheader = line.split()
                    bedstop = int(newheader[2][6:])
                    bedstart = bedstop - 1
                    # chrom should never change
                    # but just in case
                    # chrom = newheader[1][6:]
                else:
                    line = float(line)
                    print(*[chrom,bedstart,bedstop,line],sep="\t",file=out)
                    # assumes fixedstep and strepsize=1
                    bedstop += 1
                    bedstart += 1

phylopwig2bed(sys.argv[1], sys.argv[2])

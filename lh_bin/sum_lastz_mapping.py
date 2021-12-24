#!/usr/bin/env python



import argparse
import textwrap
from collections import defaultdict

def main():
    prog_name = "Another python program"
    usage = "Another python program"

    parser = argparse.ArgumentParser(
        prog=prog_name, 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("qry1", help="qry1 file")
    args = parser.parse_args()  

    qry1_file = args.qry1

#name1	zstart1	end1	name2	strand2	zstart2+	end2+	identity	idPct	coverage	covPct	cigarx-
#tig00000136	0	6386	tig00004993	+	2658608	2664786	6091/6134	99.3%	6386/60556	10.5%	171=i42=10d91=x558=10i891=x241=x44=x70=x13=x295=x65=x4=x2=x61=d15=x57=i23=i141=x59=7i622=d244=x570=x24=x104=x7=5i3=x4=x27=x6=16i87=d92=x15=x13=d18=x46=x15=x27=x=x47=x70=x147=2i177=x=i7=x110=x24=x6=x129=x10=238dx13=x56=x116=x2=x15=x220=x144=x27=x2=
#tig00000136	6342	34374	tig00004993	+	2665355	2693409	27870/27999	99.5%	28032/60556	46.3%	6=x=2x2=2x=x2=4x2=3x3=x2=x2=2x2=x5i2=x142=x135=x42=x70=x2=x=x119=x40=x14=x57=x30=x3=x115=x11=i172=x154=x22=x5=x63=x11=x73=x35=x46=x10=x4=x44=x179=x40=x74=x91=2x7=x77=x26=x48=d45=x7=x10=x17=x18=x96=2x302=x8=x24=x30=x43=x84=x27=x28=x139=x58=x18=x15=2x40=x4=x67=x33=x111=x27=x15=x37=x30=x11=19i130=x185=x=x77=x152=x147=x12=x10=x70=x91=x21=x202=x148=x219=x43=16d1131=3x870=x224=5i326=d27=x3=x12=d74=2d77=x5=i23=3d13=15i103=x27=x105=x52=x50=x20=x222=x1018=x28=8i557=x1061=4d232=x908=x1032=x358=d813=x19=x131=x1918=x263=x89=4d735=x=i5542=x897=x930=x101=x12=x1703=x375=x767=x121=x56=
    pair_dict = defaultdict(float)
    identity_cutoff = 0.97
    with open(qry1_file) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            mylist = line.rstrip().split()
            key = mylist[0] + "\t" + mylist[3]
            idenpct = float(mylist[8][:-1]) / 100
            covpct  = float(mylist[10][:-1]) / 100
            if(idenpct >= identity_cutoff):
                pair_dict[key] += covpct
    for k in pair_dict:
        print("{}\t{:.1%}".format(k, pair_dict[k]))

if __name__ == "__main__":
    main()

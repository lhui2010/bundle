#!/usr/bin/env python

import argparse
import textwrap
from collections import defaultdict
#from sympy import Interval, Union

#https://stackoverflow.com/questions/15273693/union-of-multiple-ranges
def union(data):
    """ Union of a list of intervals e.g. [(1,2),(3,4)] """
    b = []
    for begin,end in sorted(data):
        if b and b[-1][1] >= begin - 1:
            b[-1][1] = max(b[-1][1], end)
        else:
            b.append([begin, end])
    size = 0
    for begin,end in b:
        interval = end - begin + 1
        size += interval
#    print(b)
#    print(str(end) + "\t" + str(begin))
    return size
#    intervals = [Interval(begin, end) for (begin, end) in data]
#    u = Union(*intervals)
#    return [list(u.args[:2])] if isinstance(u, Interval) \
#       else list(u.args)

class Pair():
    def __init__(self, name_qry = "", name_ref = "", len_qry = "",
        len_ref = "", match_len = "", cov_qry = "", cov_ref = "", 
        loci_qry = [], loci_ref = []):
        """
        A class to hand alignment pairs
        """
        self.name_qry  = name_qry  
        self.name_ref  = name_ref  
        self.loci_qry  = [loci_qry]
        self.loci_ref  = [loci_ref]
        self.len_qry   = int(len_qry  ) 
        self.len_ref   = int(len_ref  ) 
        self.match_len = int(match_len) 
        self.cov_qry   = int(cov_qry  ) 
        self.cov_ref   = int(cov_ref  ) 

class PairDict:
    def __init__(self):
        """
        A class to hand alignment pairs
        """
        self.db = defaultdict(Pair)

    def append(self,paf_line):
        """
        Read paf file by line and store in pair file
        """
        mylist = paf_line.rstrip().split()
        name_qry = mylist[0]
        len_qry = int(mylist[1])
        name_ref = mylist[5]
        len_ref = int(mylist[6])
        if(len_qry > len_ref or name_ref == name_qry):
            #print("ERROR")
            return 0            
        else:
            cov_qry = int(mylist[3]) - int(mylist[2])
            cov_ref = int(mylist[8]) - int(mylist[7])
            match_len = mylist[9] 

            key = name_qry + "\t" + name_ref
            loci_qry = [int(mylist[2]), int(mylist[3])]
            loci_ref = [int(mylist[7]), int(mylist[8])]

            if(key not in self.db):
                self.db[key] = Pair(name_qry, name_ref, len_qry, len_ref, 
                    match_len, cov_qry, cov_ref, loci_qry, loci_ref)
            else:
                self.db[key].match_len += int(match_len)
                self.db[key].cov_qry   += int(cov_qry)
                self.db[key].cov_ref   += int(cov_ref)
                self.db[key].loci_qry.append(loci_qry)
                self.db[key].loci_ref.append(loci_ref)
                self.db[key].len_qry = int(len_qry)
                self.db[key].len_ref = int(len_ref)

    def read_paf(self,paf_file):
        """
        Read paf file by line and store in pair file
        """
# tig00000360|arrow_np1212:1-2504178      2504178 2304686 2345565 -       tig00000967|arrow_np1212:1-63923        63923   5       40704   38309   40999   0       tp:A:S  cm:i:3762       s1:i:38220      dv:f:0.0047     rl:i:12835
# tig00000360|arrow_np1212:1-2504178      2504178 2278030 2301248 -       tig00000967|arrow_np1212:1-63923        63923   40705   63912   22514   23232   0       tp:A:S  cm:i:2190       s1:i:22508      dv:f:0.0022     rl:i:12835
        with open(paf_file) as fh:
            for paf_line in fh:
                self.append(paf_line)

    def print_out(self):
        """
        print out
        """
        print("\t".join(["QryName", "RefName", "Coverage", "Identity", 
            "match_len", "cov_qry(redundant)", "cov_ref(redundant", 
            "cov_qry(after union)", "cov_ref(after union)"]))
        for k in self.db:
            line = k + "\t"
            qry_union_cov = union(self.db[k].loci_qry)
            ref_union_cov = union(self.db[k].loci_ref)
            cov_ratio_qry = qry_union_cov / self.db[k].len_qry
            identity_qry  = self.db[k].match_len / self.db[k].cov_qry

            line += "\t".join([str(cov_ratio_qry), str(identity_qry)]) + "\t"
#Error len_ref is not reported
            line += "\t".join([str(self.db[k].len_qry), str(self.db[k].len_ref), 
                str(self.db[k].match_len), str(self.db[k].cov_qry), str(self.db[k].cov_ref), 
                str(qry_union_cov), str(ref_union_cov)])
            print(line)


def main():
    prog_name = "Another python program"
    usage = "Another python program"

    parser = argparse.ArgumentParser(
        prog=prog_name, 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("qry1", help="qry1 file")
#    parser.add_argument("-f", "--flanking", default=10000, type=int, help="flanking distance default (1000)")
    args = parser.parse_args()  

    qry1_file = args.qry1
    new_pair = PairDict()
    new_pair.read_paf(qry1_file)
    new_pair.print_out()
#    flanking_distance = args.flanking




if __name__ == "__main__":
    main()

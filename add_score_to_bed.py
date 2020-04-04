import sys
from merge_bed_to_bedpe import BedIO

#==> A188B73.genome.mod.EDTA.TElib.fa.id.sort <==
#TE_00000011    DNA/DTA
#TE_00000033    DNA/DTA
#TE_00000081    DNA/DTA
#TE_00000112    DNA/DTA
#TE_00000120    DNA/DTA
#TE_00000158    DNA/DTA
#TE_00000159    DNA/DTA
#TE_00000170    DNA/DTA
#TE_00000186    DNA/DTA
#TE_00000225    DNA/DTA
#
#==> A188_gene_repeat.bed <==
#1    9    94    TE_00013772_LTR    0    +
#1    103    273    TE_00013772_LTR    0    +
#1    286    454    TE_00013772_LTR    0    +
#1    466    633    TE_00013772_LTR    0    +
#1    645    813    TE_00013772_LTR    0    +
#1    825    992    TE_00013772_LTR    0    +
#1    1004    1172    TE_00013772_LTR    0    +
#1    1184    1352    TE_00013772_LTR    0    +
#1    1364    1531    TE_00013772_LTR    0    +
#1    1543    1710    TE_00013772_LTR    0    +

def line_to_dict(input_file):
    TE_fam = {}
    with open(input_file) as fh:
        for line in fh:
            mylist = line.rstrip().split()
            TE_fam[mylist[0]] = mylist[1]
    return TE_fam

#Output
#Add family and score
#Score rule:
class teBedIO(BedIO):
    def __init__(self, filebed):
        """ init parent's attribute """
        self.bed_line = {}
        self.strand = {}
        self.score = {}
        self.chr_start_end = {}
        self.name = {}
        """ 
        Read filebed into a bed dict
        I used chromsome + start + end + strand +name +strand as key, in case sometimes one key occur in a bed multiple times
        """
        
        bed_key = ""
        with open(filebed) as fh:
            for line in fh:
                mylist = line.rstrip('\n').split('\t')
                (self.chr, self.start, self.end, gene_name, self.this_score, self.this_strand) = mylist[0:6]
                """Ignore lines with empty name"""
                if gene_name == "":
                    continue
                bed_key = self.chr + self.start + self.end + gene_name + self.this_strand
                self.bed_line[bed_key] = "\t".join(mylist[0:6])
                self.strand[bed_key] = self.this_strand
                self.score[bed_key] = self.this_score
                self.chr_start_end[bed_key] = "\t".join(mylist[0:3])
                self.name[bed_key] = gene_name
    
    def rename_with_dict(self, te_fam_dict):
        """ change TE123 to TE123#LTR/Gypsy"""    
        for k in self.name:
            name = self.name[k]
            if name in te_fam_dict:
                new_name = name + "#" + te_fam_dict[name]
                new_line = "\t".join([self.chr_start_end[k], new_name, self.score[k], self.strand[k]]) 
                self.name[k] = new_name
                self.bed_line[k] = new_line

    def change_score(self, score_dict):
        for k in self.score:
            name = self.name[k]
#Gene's score
            this_score = "1"
            if('#' in name):
                fam = name.split('#')[1]
                this_score = score_dict[fam]
            self.score[k] = this_score
            self.bed_line[k] = "\t".join([self.chr_start_end[k],self.name[k], this_score, self.strand[k]]) 

    def print_bed(self):
        for k in self.bed_line:
            print(self.bed_line[k])

def main():
    fam_score = {
    "DNA/DTA"      : "11",
    "DNA/DTC"      : "12",
    "DNA/DTH"      : "13",
    "DNA/DTM"      : "14",
    "DNA/DTT"      : "15",
    "DNA/Helitron" : "16",
    "MITE/DTA"     : "21",
    "MITE/DTC"     : "22",
    "MITE/DTH"     : "23",
    "MITE/DTM"     : "24",
    "MITE/DTT"     : "25",
    "LTR/Copia"    : "31",
    "LTR/Gypsy"    : "32", 
    "LTR/unknown"  : "33" }
#1  9   94  TE_00013772_LTR#LTR/Gypsy 31 +
    TE_fam = line_to_dict(sys.argv[1])

    beds = teBedIO(sys.argv[2])
    beds.rename_with_dict(TE_fam)
    beds.change_score(fam_score)

    beds.print_bed()

if __name__ == '__main__':
    main()

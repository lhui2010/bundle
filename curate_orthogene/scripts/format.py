from optparse import OptionParser  
from merge_bed_to_bedpe import BedIO
from Bio import SeqIO
import re


sep = "_"

def rename_fasta(prefix, file_name):
    gene_dict = SeqIO.to_dict(SeqIO.parse(file_name, "fasta"))
    print_buff = ""
    for k in gene_dict:
        new_name = re.sub(r'^', prefix + "_", k)
        print_buff += ">" + new_name + "\n" + gene_dict[k].seq.__str__() + "\n"
    return print_buff

def rename_bed(prefix, file_name):
    this_bed = BedIO(file_name)
    this_bed.rename(r'^', prefix + "_")
    print_buff = this_bed.print()
    return print_buff

def main():
    usage = """Rename gene ID with prefix
Usage: 
    {} fasta A188 A188.pep >A188.rename.pep
    {} bed A188 A188.bed > A188.rename.bed
    """.format(__file__, __file__)
    parser = OptionParser(usage=usage)
      
    (options, args) = parser.parse_args()  
#args = parser.parse_args()  

#print("options")
#print(options)
    if len(args) < 3:
        print(parser.usage)
        exit()

    if args[0] == "fasta":
        print(rename_fasta(args[1], args[2]), end='')
    elif args[0] == "bed":
        print(rename_bed(args[1], args[2]), end='')

if __name__ == "__main__":
    main()


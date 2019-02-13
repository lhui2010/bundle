import argparse
import re
import sys

#Input:
#new.genome
#old.genome
#old.gff

#Step1. extract gene_body sequences
#perl gff_genome_to_genes.pl old.gff old.genome >old.gene

#Step2. blat gene_body against new.genome
#genblast2cds_bed

#########################################################
#Input
#old.gff
#align.bed


parser = argparse.ArgumentParser(description='Liftover from aligned.bed and old gff')
parser.add_argument('GFF', type=str, nargs = 1,
                    help='old gff that need to be lifted')
parser.add_argument('BED', type=str, nargs = 1,
                    help='aligned.bed (bed of GENE(GFF\'s GENE) sequences aligned to new genome')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
args = parser.parse_args()

#print("Argument values:")
#print(args.GFF)
#print(args.BED)

def magic_loci(loci, old_start, old_end, old_strand, align_start, align_end, align_strand):
    if(align_strand == old_strand):
        transformed_loci = loci - old_start + align_start 
    else:
        transformed_loci = new_end - (loci - old_start)
    return transformed_loci 

def magic_strand(old_strand, align_strand):
    return align_strand


#Step4. read old.gff and get 1-based gene feature table
#original_scaffold_1037  maker   gene    291853  292879  .       -       .       ID=A188G39283;Name=A188G39283;Alias=maker-original_scaffold_1037-snap-gene-3.14;
#original_scaffold_1037  maker   mRNA    291853  292879  .       -       .       ID=A188G39283-t1;Parent=A188G39283;Name=A188G39283-t1;Alias=maker-original_scaffold_1037-snap-gene-3.14-mRNA-1;_AED=0.07;_QI=0|0|0.5
#original_scaffold_1037  maker   exon    291853  292724  .       -       .       ID=A188G39283-t1:2;Parent=A188G39283-t1;
#original_scaffold_1037  maker   exon    292844  292879  .       -       .       ID=A188G39283-t1:1;Parent=A188G39283-t1;
#original_scaffold_1037  maker   CDS     292844  292879  .       -       0       ID=A188G39283-t1:cds;Parent=A188G39283-t1;
#original_scaffold_1037  maker   CDS     292536  292724  .       -       0       ID=A188G39283-t1:cds;Parent=A188G39283-t1;
#original_scaffold_1037  maker   three_prime_UTR 291853  292535  .       -       .       ID=A188G39283-t1:three_prime_utr;Parent=A188G39283-t1;

old_start = {}
old_strand = {}

with open (args.GFF[0]) as fh_gff:
    for line in fh_gff:
        (contig, tmp, feature_type, 
            start, end, tmp2, strand, 
            phase, feats) = line.split()
        if(feature_type == "gene"):
            gene = re.match("ID=(.*?);", feats)[1]
            old_start[gene] = int(start)
#            old_end[gene] = end
            old_strand[gene] = strand

#Step3. get new_gene_offset
#Assuming bed is aligned using oriented (5'->3') gene sequence
#def main():
align_start = {}
align_end = {}
align_strand = {}
align_contig = {}
offset = {}
count_bad =0
count_all =0
with open (args.BED[0]) as fh_bed:
    for line in fh_bed:
        (contig, start, end, gene, tmp, strand) = line.split()
        align_start[gene] = int(start)
        align_end[gene] = int(end)
        align_strand[gene] = strand
#strand on new coordinate system is the same with aligned direction
        align_contig[gene] = contig
        if(align_strand[gene] == old_strand[gene]):
            offset[gene] = align_start[gene] - old_start[gene]
        else:
            offset[gene] = 0- align_end[gene] - old_start[gene]

#Step5. read old.gff again and output new.gff
with open (args.GFF[0]) as fh_gff:
    for line in fh_gff:
        (contig, tmp, feature_type, 
            start, end, tmp2, strand, 
            phase, feats) = line.split()
#magic formular
        if(feature_type == "gene"):
            gene = re.match(r'ID=(.*?);', feats)[1]
        else:
            gene = re.match(r'.*Parent=(.*?);', feats)[1]
        gene = re.sub('-.*','',gene)
        count_all +=1
        if(not offset.__contains__(gene)):
            count_bad +=1
            sys.stderr.write(gene+":\t" +line)
            continue
        new_start = abs(offset[gene] + int(start))
        new_end = abs(offset[gene] + int(end))
        new_strand = align_strand[gene]
        new_chr = align_contig[gene]
        
        new_line = "\t".join([new_chr, ".", feature_type,
            str(new_start), str(new_end), ".", new_strand,phase, feats.rstrip()])
        print(new_line)

sys.stderr.write("count_all:\t" + str(count_all)+"\n")
sys.stderr.write("count_bad:\t" + str(count_bad)+"\n")

#Output:
#new.gff

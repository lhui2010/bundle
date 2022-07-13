from Bio import AlignIO
alignment = AlignIO.read(open("MtSHR_ortho.txt.fa.aln.rename"), "fasta")
print(format(alignment, "maf"))

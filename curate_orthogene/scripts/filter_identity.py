#A188G00063-t1	A188G00063-t1-R1-1-A1	1.0	A188G00063-t1	B73_Zm00001d033554_T001	0.2900763358778626
#A188G00074-t1	A188G00074-t1-R1-1-A1	1.0	A188G00074-t1	B73_Zm00001d033544_T001	0.9891745602165088
#A188G11692-t1	A188G11692-t1-R3-3-A1	0.9787234042553191	A188G11692-t1	B73_Zm00001d049360_T001	0.9432624113475178
#A188G24714-t1	A188G24714-t1-R2-2-A1	0.8160919540229885	A188G24714-t1	B73_Zm00001d052930_T002	0.4482758620689655
import sys

identity_cutoff = 0.3

if(len(sys.argv) > 2):
    identity_cutoff = float(sys.argv.pop(1))

with open(sys.argv[1]) as fh:
    for line in fh:
        (REF_gene, ortho_genblast, iden_genblast, 
        REF_gene2, ortho_geneset, iden_geneset) = line.rstrip().split()
        target_ortho = ortho_geneset
        if(float(iden_genblast) - float(iden_geneset) > identity_cutoff):
            target_ortho =  ortho_genblast
            #print(line.rstrip())
            print("{}\t{}".format(REF_gene, target_ortho))
        

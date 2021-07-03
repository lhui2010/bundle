## Gene annotation with MAKER

### RepeatMasker

* Input
    - Genome
* Output
    - Repeat.gff (GFF format of repeat)
    - Repeat.tbl (Summary of repeat ratio)
    - Softmasked genome (complex repeat only, for maker and braker annotation)

```
# Step1. Run repeatmasker with both homology and de novo approach
python -m iga.annotation.repeat repeatmasker --species Viridiplantae elumb.contig.fa --threads 10
## Output:
## workdir_repeatmask_elumb.contig.fa/species_lib.out/elumb.contig.fa.cat.gz  
## workdir_repeatmask_elumb.contig.fa/species_lib.out/elumb.contig.fa.out
## workdir_repeatmask_elumb.contig.fa/species_lib.out/elumb.contig.fa.masked  
## workdir_repeatmask_elumb.contig.fa/species_lib.out/elumb.contig.fa.tbl

## The following step is slow  
python -m iga.annotation.repeat repeatmasker --denovo T elumb.contig.fa --threads 10
## Output:
## workdir_repeatmask_elumb.contig.fa/custom_lib.out/elumb.contig.fa.cat.gz  
## workdir_repeatmask_elumb.contig.fa/custom_lib.out/elumb.contig.fa.out
## workdir_repeatmask_elumb.contig.fa/custom_lib.out/elumb.contig.fa.masked  
## workdir_repeatmask_elumb.contig.fa/custom_lib.out/elumb.contig.fa.tbl

# Step2. Combine the two result 

python -m iga.annotation.repeat post_repeatmasker workdir_repeatmask_elumb.contig.fa --genome elumb.contig.fa

```

### Prepare RNA-seq evidence

#### Get flnc reads (fasta format)

```
for i in *flcdna*bam; do python -m iga.annotation.isoseq isoseq_pb $i elumb.flcdna.pb.primer.fasta; done
```

#### Assemble RNA-Seq reads

```
python -m iga.annotation.rnaseq reads_align_assembly "zenin.ssRNA.ZnY5_1.clean.fq.gz zenin.ssRNA.ZnY5_2.clean.fq.gz" zeins.contig.fa --threads 2
```

#### RNA-Seq fasta to gff

```
# daemon will not exist until the job finished
python -m iga.annotation.maker fastq2gff est.fasta genome.fasta
```

### Prepare Protein evidence

```
python -m iga.annotation.maker prep_genblast elumb.contig.fa.masked.fa pep/elumb.homologous_pep.fa
```

* Input: 
    - rna_fasta.bam
    - protein_evidence.fasta
    - Genome (soft masked of complex region)
* OutputHMM: test_datisca/species/Sp_5/


### BRAKER


```
wd=test_datisca

if [ -d $wd ]; then
    rm -r $wd
fi

REF=Datisca_v2.masked.fa
PEP=total_ortho.pep
ESTBAM=Trinity-GG.fasta.bam

( time braker.pl --genome=${REF} --prot_seq=${PEP} --prg=gth --bam=${ESTBAM} --gth2traingenes \
--softmasking --workingdir=$wd ) &> $wd.log
```

### MAKER

```

REF=elumb.contig.fa.masked.fa
FLNCESTGFF=flnc_rna.gff
ESTGFF=total_est.gff
CDNAFASTA=flnc_rna.fasta
PEPGFF=pep.gff
REPEATGFF=repeat.gff


#-+-First Round
ROUND=1
#python -m iga.annotation.maker deploy_augustus
python -m iga.annotation.maker maker_run         ${REF} ${FLNCESTGFF} ${PEPGFF} ${REPEATGFF} --cpus 2
python -m iga.annotation.maker maker_check_resub ${REF}_R1
python -m iga.annotation.maker maker_collect     ${REF}_R1
python -m iga.annotation.maker maker_train       ${REF}_R1 --cdna_fasta ${CDNAFASTA}  --snap 'T' --augustus F
cd ${REF}_R${ROUND}
python -m iga.assembly.assess busco --mode prot total.all.maker.proteins.fasta
cd ..

#-+-Second Round
PREV=1
ROUND=2
# python -m iga.annotation.maker deploy_augustus
python -m iga.annotation.maker maker_run         ${REF} ${FLNCESTGFF} ${PEPGFF} ${REPEATGFF} --round ${ROUND} --augustus_species ${REF}_R${PREV}_direct --snap_h
mm ${REF}_R${PREV} --update "alt_splice=1"
python -m iga.annotation.maker maker_check_resub ${REF}_R${ROUND}
python -m iga.annotation.maker maker_collect     ${REF}_R${ROUND}
python -m iga.annotation.maker maker_train       ${REF}_R${ROUND} --cdna_fasta ${CDNAFASTA} --augustus F
cd ${REF}_R${ROUND}
python -m iga.assembly.assess busco --mode prot total.all.maker.proteins.fasta
cd ..

#-+-Third Round #augustus is directly trained
ROUND=3
PREV=2
# python -m iga.annotation.maker deploy_augustus
python -m iga.annotation.maker maker_run         ${REF} ${ESTGFF} ${PEPGFF} ${REPEATGFF} --round ${ROUND} --augustus_species ${REF}_R${PREV}_direct --snap_hmm $
{REF}_R${PREV} --update "trna=1;alt_splice=1"
#--queue Q64C1T_X4
python -m iga.annotation.maker maker_check_resub ${REF}_R${ROUND}
python -m iga.annotation.maker maker_collect     ${REF}_R${ROUND}
cd ${REF}_R${ROUND}
python -m iga.assembly.assess busco --mode prot total.all.maker.proteins.fasta
cd ..


python -m iga.annotation.maker pasa_refine ref.fa ../flnc_rna.fasta genome.maker.gff --use_grid T

$bsub512 "python -m iga.annotation.maker maker_rename_gff pasa_raw.gff3"
grep trna ../genome.maker.gff > trna.gff

```
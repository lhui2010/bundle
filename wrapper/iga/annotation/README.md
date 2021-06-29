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
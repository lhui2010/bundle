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

species=Viridiplantae
cd {0}
mkdir -p Full_mask
## unzip
gunzip *lib.out/*.cat.gz
cat *lib.out/*.cat >full_mask.cat
## to mask.out
ProcessRepeats -species $species full_mask.cat
## create GFF3
rmOutToGFF3.pl full_mask.out > full_mask.gff3

## isolate complex repeats
grep -v -e "Satellite" -e ")n" -e "-rich" full_mask.gff3 \
     > full_mask.complex.gff3

## reformat to work with MAKER
cat full_mask.complex.gff3 | \
 perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' \
         > full_mask.complex.reformat.gff3

echo "Repeat GFF file is located in "
echo "$PWD/full_mask.complex.reformat.gff3"
```

# Step3. Get soft masked genome
 
```
bedtools maskfasta -soft -fi elumb.contig.fa -bed full_mask.complex.reformat.gff3 \
-fo elumb.contig.masked.fa
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
python -m iga.annotation.maker fastq2gff 
```

### Prepare Protein evidence

* Input: 
    - rna_fasta.bam
    - protein_evidence.fasta
    - Genome (soft masked of complex region)
* OutputHMM: test_datisca/species/Sp_5/

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


### BRAKER


### MAKER
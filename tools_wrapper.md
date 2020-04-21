
#### conda
export PS1="(base) \[\033]2;\h:\u $PWD\007\033[33;1m\]\u@\h \033[35;1m\t\n\033[0m\[\033[36;1m\]$PWD\[\033[0m\]\n\[\e[32;1m\]$\[\033[0m\]"
source ~/lh/anaconda2/etc/profile.d/conda.sh
conda activate jcvi


#### last align cds
lastdb sorghum sorghum.cds
lastal -u 0 -P 56 -i3G -f BlastTab sorghum A188.cds >A188.sorghum.last

#### RAxML
for i in func_gene.fna
do
    qsub -V -b y -N output -cwd "raxmlHPC-PTHREADS-AVX2 -f a -x  32431 -p 43213 -# 100 -m GTRGAMMA -T 40 -s ${i} -n ${i%.*} "
    qsub -V -b y -N output -cwd "raxmlHPC-PTHREADS-AVX2 -f a -x  32431 -p 43213 -# 100 -m PROTGAMMAWAGF -T 40 -s ${i} -n ${i%.*} "
    qsub -V -b y -N output -cwd ' raxml-ng --all --msa total.c2h2.YJ.pep.aln.fa --model LG+G8+F --prefix YJpep --msa-format FASTA'
done

#### Blast
makeblastdb -in input.fa -dbtype nucl   prot
blastp -subject sub.fa -query qry.fa -out out.bln -evalue 1e-5 -outfmt 6 -num_threads 16

#### genewise 
```
#Note the prediction do not include the last stop codon
genewise -cdna -pseudo Zm00001d042922_T002.pep target.fa >Zm00001d042922_T002.fa.vs.target.gene_wise
```
#### bwa
bwa index consensus.fasta
bwa mem consensus.fasta -t 40 read1.gz read2.gz >bwa.sam 2>bwa.err
bwa mem consensus.fasta -t 40 read1.gz >bwa.sam 2>bwa.err

#### samtools
samtools view -hbS xx.sam >ss.bam
samtools mpileup -uvf M445.chr.fa M441_sequences_to_M445_genome.bam  >tmp.vcf
samtools tview -d T -p chr01:2473903 subset.bam

### bsub

bswitch target_quename 5309(job_id)
bsub -m "hostA hostD hostB" myjob
https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.3/lsf_admin/job_view_all.html
bjobs -u all
bqueues -u yitingshuang
bqueues -l Q104C512G_X4 |grep HOSTS
bqueues -l Q64C1T_X4 |grep HOSTS
bqueues -l Q48C2T_X1 |grep HOSTS
bqueues -l Q64C3T_X1 |grep HOSTS

#### Monitor jobs
while true; do date; lsload |grep "HOST\|node02\|node03\|node04\|node05"; sleep 1m; done

### fastqc
$bsub512 -J QC  'mkdir fastqc_dir && perl  /ds3200_1/proc/FastQC/fastqc -t 60 *gz -o fastqc_dir'

### snpEFF

GENOME=/lustre/home/liuhui/project/A188/maize_ortho/A188.genome
GFF=/lustre/home/liuhui/project/A188/maize_ortho/A188.gff
DBNAME=A188

set -euxo pipefail

mkdir -p ${SE_HOME}/data/${DBNAME}
cp ${GFF} ${SE_HOME}/data/${DBNAME}/genes.gff
cp ${GENOME} ${SE_HOME}/data/${DBNAME}/sequences.fa
cd ${SE_HOME}
echo ${DBNAME}.genome : ${DBNAME} >>snpEff.config
java -jar ${SE_HOME}/snpEff.jar build -gff3 -v ${DBNAME}

java -jar ${SE_HOME}/snpEff.jar eff Gossypium_arborium /ds3200_1/users_root/yitingshuang/lh/projects/polyploid_evol/03.GO_enrichment/Ga09G1341.pep.vcf >snpEff_genes.txt.info

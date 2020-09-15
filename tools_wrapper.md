
#### conda
export PS1="(base) \[\033]2;\h:\u $PWD\007\033[33;1m\]\u@\h \033[35;1m\t\n\033[0m\[\033[36;1m\]$PWD\[\033[0m\]\n\[\e[32;1m\]$\[\033[0m\]"
source ~/lh/anaconda2/etc/profile.d/conda.sh
conda activate jcvi

#### conda env
$conda env list
# conda environments:
#
                         /ds3200_1/users_root/yitingshuang/lh/anaconda2
                         /ds3200_1/users_root/yitingshuang/lh/anaconda2/envs/jcvi
                         /ds3200_1/users_root/yitingshuang/lh/anaconda2/envs/jitterbug
                         /ds3200_1/users_root/yitingshuang/lh/anaconda2/envs/longshot
base                     /ds3200_1/users_root/yitingshuang/lh/anaconda3
EDTA                     /ds3200_1/users_root/yitingshuang/lh/anaconda3/envs/EDTA
busco                    /ds3200_1/users_root/yitingshuang/lh/anaconda3/envs/busco
deepvariant              /ds3200_1/users_root/yitingshuang/lh/anaconda3/envs/deepvariant
falcon                   /ds3200_1/users_root/yitingshuang/lh/anaconda3/envs/falcon
hyphy                    /ds3200_1/users_root/yitingshuang/lh/anaconda3/envs/hyphy
longshot                 /ds3200_1/users_root/yitingshuang/lh/anaconda3/envs/longshot
r_env                 *  /ds3200_1/users_root/yitingshuang/lh/anaconda3/envs/r_env (Used for smudge_plot)
syri                     /ds3200_1/users_root/yitingshuang/lh/anaconda3/envs/syri
trinity                  /ds3200_1/users_root/yitingshuang/lh/anaconda3/envs/trinity



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

#### maker mapid
cp assemblyv1.est2genome.gff assemblyv1.est2genome.gff.bak
maker_map_ids --prefix "corv1" --suffix '-t' --iterate 1 assemblyv1.est2genome.gff >corv1.id.map
map_gff_ids corv1.id.map  assemblyv1.est2genome.gff

#### genewise 
```
#Note the prediction do not include the last stop codon
genewise -cdna -pseudo Zm00001d042922_T002.pep target.fa >Zm00001d042922_T002.fa.vs.target.gene_wise
```
#### bwa
REF=xx
FASTQ1=xx
FASTQ2=xx
THREADS=80
bwa index $REF
bwa mem -t ${THREADS} $REF ${FASTQ1} ${FASTQ2} \
| samtools view -@ ${THREADS} -bS - | samtools sort -@ 80 - > ${REF}.bam
samtools index input/${REF}.bam
samtools flagstat -@ ${THREADS} ${REF}.bam >${REF}.bam.stat

bwa index consensus.fasta
bwa mem consensus.fasta -t 40 read1.gz read2.gz >bwa.sam 2>bwa.err
bwa mem consensus.fasta -t 40 read1.gz >bwa.sam 2>bwa.err

#### samtools
samtools view -hbS xx.sam >ss.bam
samtools mpileup -uvf M445.chr.fa M441_sequences_to_M445_genome.bam  >tmp.vcf
samtools tview -d T -p chr01:2473903 subset.bam

#### Get unmapped
samtools view -f 4 falcon_v340.fasta.bam |head -4000 |sam2fq.pl  |fq2fa.pl - >unmapped.fa

### bsub

bswitch target_quename 5309(job_id)
bsub -m "hostA hostD hostB" myjob
https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.3/lsf_admin/job_view_all.html
bjobs -u all
bqueues -u yitingshuang
bqueues -l Q104C512G_X4 |grep HOSTS
    HOSTS:  node10 node11 node12 node13
bqueues -l Q64C1T_X4 |grep HOSTS
    HOSTS:  node02 node03 node04 node05
bqueues -l Q48C2T_X1 |grep HOSTS
    HOSTS:  node07
bqueues -l Q64C3T_X1 |grep HOSTS
    HOSTS:  node06

#### unavalable queues
bqueues -l Q104C512G_X2 |grep HOSTS
    HOSTS:  node08 node09
bqueues -l Q88C6T_X1 |grep HOSTS
    HOSTS:  node01
bqueues -l Q224C12T_X1 |grep HOSTS
    HOSTS:  node15
bqueues -l QGPU |grep HOSTS
    HOSTS:  node14



#### Monitor jobs
while true; do date; lsload |grep "HOST\|node02\|node03\|node04\|node05"; sleep 1m; done
while [ `bjobs 127421 2>/dev/null |head -1 |awk '{print $1}'` == JOBID ]; do sleep 1m; done

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

### minimap2

#### Error correction
minimap2 -t${THREADS} -ax map-pb -r2k ${REF} ${READS} | samtools sort -@${THREADS} >subreads_on_mt.bam

#### ALign isoseq

minimap2 -t 30 -ax splice -uf --secondary=no -C5 -O6,24 -B4 \
   hg38.fasta hq_isoforms.fasta \
   > hq_isoforms.fasta.sam \
   2> hq_isoforms.fasta.sam.log


### R devtools
devtools::install_deps(dependencies = TRUE)

### dotPlot
minimap2 -t 10 -x asm5 ${REF} ${contig} > ${contig}.paf
/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/nextdenovo/dotPlotly/pafCoordsDotPlotly.R  -i ${contig}.paf -o ${contig}  -l -p 6 -k 12


### lastz
lastz --gfextend --chain --gapped --ambiguous=n --format=mapping --identity=97  --continuity=1  ${REF} ${QRY}

### trf
trf yoursequence.txt 2 7 7 80 10 50 500 -f -d -m 

### canu

ROOT=$PWD
CANU=/ds3200_1/users_root/yitingshuang/lh/bin/canu/canu-2.0/Linux-amd64/bin/canu
PREFIX=Coriaria_unmapped
WORKDIR="workdir_"${PREFIX}
INPUT=${ROOT}/correads_on_purged.bam.unmap.sam.fa
#pacbio or nanopore
TYPE=pacbio-hifi

mkdir -p ${WORKDIR}
#polyploid_param
# corOutCoverage=200 "batOptions=-dg 3 -db 3 -dr 1 -ca 500 -cp 50" \

${CANU} \
 -d ${WORKDIR} -p ${PREFIX}  \
 -s threads.config \
genomeSize=1m \
useGrid=true \
minReadLength=1000 \
minOverlapLength=500 \
-${TYPE} \
${INPUT} \
  >>${WORKDIR}/assemble.log 2>>${WORKDIR}/assemble.err

#### canu threads.config
more threads.config 

maxMemory=500g
maxThreads=64

merylThreads=64
merylMemory=200

cormhapThreads=8
cormhapConcurrency=1
obtmhapThreads=8
utgmhapThreads=8

cnsThreads=8
corThreads=8

corovlThreads=8
obtovlThreads=8
utgovlThreads=8

### docker
sudo yum install -y yum-utils
sudo yum-config-manager \
    --add-repo \
        https://download.docker.com/linux/centos/docker-ce.repo
sudo yum install docker-ce docker-ce-cli containerd.io
sudo echo "{
      "registry-mirrors": ["https://docker.mirrors.ustc.edu.cn/"]
}" >> /etc/docker/daemon.json
sudo systemctl start docker
sudo docker run hello-world

### apollo
docker pull gmod/apollo
docker run -d -it --privileged --rm -p 9999:8080 -v /tmp/apollo_data gmod/apollo

#### To visit
ServerIP:9999
user: admin@local.host
pw: password
#### New organism
faToTwoBit genome.fa genome.2bit


### LTR_retriever
gt suffixerator -db genome.fa -indexname genome.fa -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index genome.fa -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > genome.fa.harvest.scn
LTR_FINDER_parallel -seq genome.fa -threads 10 -harvest_out -size 1000000 -time 300
cat genome.fa.harvest.scn genome.fa.finder.combine.scn > genome.fa.rawLTR.scn

#### To run LTR_retriever:
LTR_retriever -genome genome.fa -inharvest genome.fa.rawLTR.scn -threads 10 [options]

#### To run LAI:
LAI -genome genome.fa -intact genome.fa.pass.list -all genome.fa.out [options]


#### novogene
lnd -u -p
lnd list 
lnd cp -d oss://CP2019100800080/H101SC20072349/KY_nuohe_JK/X101SC20072349-Z01/X101SC20072349-Z01-J001/2.cleandata .


### fastqdump

[link](https://edwards.sdsu.edu/research/fastq-dump/)
fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR_ID

### slurm

#### submit jobs
sbatch submit.sh

head ~/Test_kmer/submit.sh 
#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64

#### view jobs
squeue

#### view  node stat
$sinfo -n p2407 -o '%c %m %O %T'

#### gff2bed
GFF=falcon_peps.gff
grep -v "^#" ${GFF} |sed 's/;.*//; s/ID=//' | awk '$3=="protein_match"'| awk '{print $1"\t"$4"\t"$5"\t"$9"\t"$6"\t"$7}' >${GFF}.bed

#### HYPHY

##### Command

HYPHYMP OG0008766.bs > OG0008766.hyphy_log"

##### FileType

==============================

$more *.bs
fileToExe = "/lustre/home/liuhui/bin/anaconda3/lib/hyphy/TemplateBatchFiles/SelectionAnalyses/aBSREL.bf";


inputRedirect = {};
inputRedirect["01"]="Universal";
inputRedirect["02"]="/lustre/home/liuhui/project/buzzo/OrthoFinder/running/maize_v1.1/Alignments/OG0008766/OG0008766.fna";
inputRedirect["03"]="/lustre/home/liuhui/project/buzzo/OrthoFinder/running/maize_v1.1/Alignments/OG0008766/OG0008766.nwk";
inputRedirect["04"]="All";

ExecuteAFile( fileToExe, inputRedirect);

==============================

$head /lustre/home/liuhui/project/buzzo/OrthoFinder/running/maize_v1.1/Alignments/OG0008766/OG0008766.fna /lustre/home/liuhui/project/buzzo/OrthoFinder/running/maize_v1.1/Alignments/OG0008766/OG0008766.nwk
==> /lustre/home/liuhui/project/buzzo/OrthoFinder/running/maize_v1.1/Alignments/OG0008766/OG0008766.fna <==
>A188G21859-t1
------------------------------------------------------------
------------------------------------------------------------
------------------------------------ATGTCCACCGCCGGGGACCCTTCC
CGCCTCTCCGGCGAGTCCTCGCCGTCGTCCTCCACGTCCTCCGGCTCCTCCTCCCACTCC
---TCTGGCGCCGCCGATGCCGCCGCCACCAACCTCGCCCTGACAGCACCAACCTCCGCC
CTCGCCGATGACACAGACGCCGATGCCCCCACCTCCCCGCGCGTGGGGACGTACTTTGAG
ACCGAGGACGACGCGTACGAGTTCTACAAGGCCTACGCGGCCCGTCTCGGCTTCGTCGTC
CGCAAGTCCAACAAGTCCAAGAACTCACGGCACACCGTCACCCGCCGCCTCTTCGTCTGC
TCCAAGCAGGGCTTCCGCCAGGAGCCCAAGAAGCCCCAGGACGAAACCGCAGGCTCCGGA

==> /lustre/home/liuhui/project/buzzo/OrthoFinder/running/maize_v1.1/Alignments/OG0008766/OG0008766.nwk <==
(Zm00008a028894_P01:0.014261724,sorghum|KXG36267:0.054149033,(Zm00001d021545_T006:0.000000005,(A188G21859-t1:0.003568854,(PWZ14309.1:0.0,Zm00004b036824_P001:0.0):0.018688784)0.944:0.011343172)0.000:0.000000005);

==============================

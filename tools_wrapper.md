
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
blastp -db Andira_inermis_Pap.pep -query Andira_inermis_Pap.pep -out Andira_inermis_Pap.pep.yangya.bln -evalue 1e10 -outfmt "6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads 10

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
samtools index ${REF}.bam
samtools flagstat -@ ${THREADS} ${REF}.bam >${REF}.bam.stat

bwa index consensus.fasta
bwa mem consensus.fasta -t 40 read1.gz read2.gz >bwa.sam 2>bwa.err
bwa mem consensus.fasta -t 40 read1.gz >bwa.sam 2>bwa.err

##### This is a short post on how to remap short reads in an aligned BAM using bwa-mem. My recommendation is (requiring bash)

samtools collate -Oun128 in.bam | samtools fastq -OT RG,BC - \
  | bwa mem -pt8 -CH <(samtools view -H in.bam|grep ^@RG) ref.fa - \
  | samtools sort -@4 -m4g -o out.bam -


#### check contamination

```
# bwa mapping

samtools view -f 4 falcon_v340.fasta.bam |head -4000 |sam2fq.pl  |fq2fa.pl - >unmapped.fa

seqtk seq -A unmapped.fq > unmapped.fa

diamond blastx --db nr.dmnd --query unmapped.fq --out unmapped.fq.tab

blastn -query unmapped.fa -out out.xml -max_target_seqs 1 -outfmt 5 -db ~/lh/database/nt -num_threads 2 -evalue 1e-5
```

#### samtools
samtools view -hbS xx.sam >ss.bam
samtools mpileup -uvf M445.chr.fa M441_sequences_to_M445_genome.bam  >tmp.vcf
samtools tview -d T -p chr01:2473903 subset.bam
samtools view ccs.bam | awk '{print ">"$1"\n"$10}' > ccs.fa

#### bcftools
bgzip Zenia_insignis.fasta.bam.vcf 
tabix Zenia_insignis.fasta.bam.vcf.gz 
cat Zenia_insignis.fasta | vcf-consensus  Zenia_insignis.fasta.bam.vcf.gz  > Zenia_insignis.fasta.bam.vcf.gz.fa 

#### Get unmapped
samtools view -f 4 falcon_v340.fasta.bam |head -4000 |sam2fq.pl  |fq2fa.pl - >unmapped.fa

### bsub

```
bswitch target_quename 5309(job_id)
bsub -m "hostA hostD hostB" myjob
bsub -J ${p:0:3}
https://www.ibm.com/support/knowledgecenter/en/SSETD4_9.1.3/lsf_admin/job_view_all.html

#Change Priority
#You can move jobs that are pending with the btop command.
btop job_ID | "job_ID[index_list]" [position]
#If you add [position] it means that the job will be put at that place in the queue.
#By default, LSF dispatches jobs in a queue in the order of their arrival (that is, first come, first served), subject to availability of suitable server hosts.


#Display Detailed Job information (including finished)
bjobs -l jobs_id
#Display All jobs
bjobs -u all
#View Job realtime STDOUT
bpeek 206654
#View all information about queues
bqueues -l
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
```

#### Monitor jobs
while true; do date; lsload |grep "HOST\|node02\|node03\|node04\|node05"; sleep 1m; done
while [ `bjobs 127421 2>/dev/null |head -1 |awk '{print $1}'` == JOBID ]; do sleep 1m; done

### sbatch
```

#Find idle nodes
sinfo -N --states=idle
```

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

#### for telomere discovery
trf CORNE_v1.0.chr.fa.ends.fa 2 9 9 80 10 100 7 -h

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

#### Best Practice

```
sbatch --get-user-env -n 1 -N 1 -c 30 -p xhacnormalb bsub.repeat_masker_Chamaecrista_fasciculata.553526.sh

```

#### submit jobs
sbatch submit.sh

head ~/Test_kmer/submit.sh 
#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64

#### interactive jobs
srun -p com -w comput[1-2] -N 2 -n 40 -t 20 A.exe
交互式提交A.exe程序。如果不关心节点和时间限制，可简写为srun -N 2 -n 40 A.exe
-p com指定提交作业到com队列；
-w comput[1-2] 指定使用节点comput[1-2]；
-N 2 指定使用2个节点；
-n 40 指定进程数为40；
-t 20 指定作业运行时间限制为20分钟。

#### ssh to node
salloc分配模式作业提交

#### view jobs
squeue
https://hpc.sicau.edu.cn/syzn/slurm.htm
squeue -j 123456
查看作业号为123456的作业信息
squeue -u paratera
查看超算账号为 paratera的作业信息
squeue –p com
查看提交到com队列的作业信息
squeue -w comput1
查看使用到comput1节点的作业信息

#### view node statistics
scontrol show node $i
scontrol show job 123456

#### view  node stat
指定显示节点comput1的使用情况
sinfo -n p2407 -o '%c %m %O %T'
指定显示队列com情况
sinfo -p com

### gff2bed
GFF=falcon_peps.gff
grep -v "^#" ${GFF} |sed 's/;.*//; s/ID=//' | awk '$3=="protein_match"'| awk '{print $1"\t"$4"\t"$5"\t"$9"\t"$6"\t"$7}' >${GFF}.bed

#### HYPHY

##### Command

```
qsub -pe smp 5 -V -b y -N $g -cwd "cd $WORKDIR/$g && t_coffee ${g}.pep -mode fmcoffee > ${g}.aln.out && cp ${g}.aln ${g}.pep.aln &&  pal2nal.pl ${g}.pep.aln ${g}.cds   >${g}.paml_aln && trimal -in ${g}.paml_aln -out ${g}.fna     -fasta && trimal -in ${g}.pep.aln -out ${g}.pep.fna -fasta && fasttree ${g}.pep.fna >${g}.nwk && touch -a ${g}.fna.ABSREL.json && mv ${g}.fna.ABSREL.json ${g}.fna.ABSREL.json.bak && HYPHYMP ${g}.bs > ${g}.hyphy_log"
```

HYPHYMP OG0008766.bs > OG0008766.hyphy_log"

grep "^* Zm00008" aBSREL.txt >aBSREL.txt.PH207

##### FileType

```

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

```


#### liftover gff, agp, and build chr.fa

```
perl -w juicer_assembly_to_ragoo_ordering.pl falcon_v340_sgs_polish.1.review.assembly
python make_agp.py orderings.fofn *fai 100 |sed 's/_RaGOO//' > CORNE_v1.0.agp
python lift_over.py  ../CORNE.contig.gff orderings.fofn CORNE.contig.fa.fai |sed 's/_RaGOO//' > CORNE.chr.gff
build_fa_from_agp.pl CORNE.contig.fa CORNE_v1.0.agp > CORNE_v1.0.chr.fa
```

#### mask fasta

`bedtools maskfasta -soft -fi ../../input/falcon_ref.fa -bed full_mask.complex.reformat.gff3 -fo falcon_ref.masked.fa`

#### genblast

merge_genblast.sh */*.gff > ../new_genblast_0921.gff
filter_genblast.py new_genblast_0921.gff > new_genblast_0921.slim.gff
sed 's/transcript/protein_match/; s/coding_exon/match_part/' new_genblast_0921.slim.gff > new_genblast_0921.slim.MAKER.gff

#### EVM

```
$EVM_HOME/EvmUtils/partition_EVM_inputs.pl --genome ${REF} --gene_predictions ${ROOT}/${GFF} \
    --transcript_alignments ${EST} \
    --segmentSize 100000 \
    --overlapSize 10000 --partition_listing partitions_list.out


    #--protein_alignments ${PEP} \
$EVM_HOME/EvmUtils/write_EVM_commands.pl --genome ${REF} --weights `pwd`/weights.txt \
    --gene_predictions $ROOT/${GFF}  \
    --transcript_alignments ${EST} \
    --output_file_name evm.out  --partitions partitions_list.out >  commands.list

$EVM_HOME/EvmUtils/execute_EVM_commands.pl commands.list | tee run.log
 #
$EVM_HOME/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
 #
$EVM_HOME/EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome ${REF}
```

#### lastz plot

```
echo "#name1	zstart1	end1	name2	strand2	zstart2+	end2+	identity	idPct	coverage	covPct	cigarx-	chr" > total.lastz.out
for i in `seq 10`; do perl -e '$id = shift @ARGV;  while(<>){next if (/^#/); chomp; $_.="\tchr$id\n";print;}' $i lastz.$i.txt >> total.lastz.out; done
lastz_ggplot.total.R total.lastz.out
```

#### bam2fastq

```
bam2fastq -o myEcoliRuns m54008_160330_053509.subreads.bam m54008_160331_235636.subreads.bam
```

#### clean fastq

```
for SAMPLE in `cat total_sample.txt`
do
    LEFT=${SAMPLE}_1.fq
    RIGHT=${SAMPLE}_2.fq
    echo "$RIGHT
$LEFT" > $SAMPLE.fastuniq
   # fastuniq -i $SAMPLE.fastuniq -o $RIGHT.dedup -p $LEFT.dedup && \
    fastx_trimmer -i $RIGHT.dedup -o ${RIGHT}.dedup.trim.gz -z -Q 33 -f 46 -l 145 && \
    fastx_trimmer -i $LEFT.dedup -o ${LEFT}.dedup.trim.gz -z -Q 33 -f 4 && \
    java -jar  /ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/trimmomatic-0.38.jar PE \
  -phred33 ${LEFT}.dedup.trim.gz ${RIGHT}.dedup.trim.gz \
   ${SAMPLE}_1.clean.fq ${SAMPLE}_1.clean.unpair.fq \
   ${SAMPLE}_2.clean.fq ${SAMPLE}_2.clean.unpair.fq \
   ILLUMINACLIP:/ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/adapters/P5P7-PE.fa:0:30:10 \
   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75 &
done

for SAMPLE in `cat total_sample.txt`
do
    LEFT=${SAMPLE}_1.fq
    RIGHT=${SAMPLE}_2.fq
    echo "$RIGHT
$LEFT" > $SAMPLE.fastuniq
   # fastuniq -i $SAMPLE.fastuniq -o $RIGHT.dedup -p $LEFT.dedup && \
    fastx_trimmer -i $RIGHT.dedup -o ${RIGHT}.dedup.trim.gz -z -Q 33 -f 46 -l 145 && \
    fastx_trimmer -i $LEFT.dedup -o ${LEFT}.dedup.trim.gz -z -Q 33 -f 4 && \
    java -jar  /ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/trimmomatic-0.38.jar PE \
  -phred33 ${LEFT}.dedup.trim.gz ${RIGHT}.dedup.trim.gz \
   ${SAMPLE}_1.clean.fq ${SAMPLE}_1.clean.unpair.fq \
   ${SAMPLE}_2.clean.fq ${SAMPLE}_2.clean.unpair.fq \
   ILLUMINACLIP:/ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/adapters/P5P7-PE.fa:0:30:10 \
   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75 &
done
```

#### eggnog

emapper.py --target_taxa "Viridiplantae" -o xx xx.pep


#### gatk

##### make alternative reference
gatk CreateSequenceDictionary -R frankia_cn3_genome.fa
gatk IndexFeatureFile -F LIU.rep3.tuxedo.bam.vcf
gatk FastaAlternateReferenceMaker -R frankia_cn3_genome.fa -V LIU.rep3.tuxedo.bam.vcf -O frankia_cn3_genome.rep3.fa

#### bcftools consensus reference

bgzip LIU.rep3.tuxedo.bam.format.vcf
tabix LIU.rep3.tuxedo.bam.format.vcf.gz
bcftools consensus -f frankia_cn3_genome.fa -o frankia_cn3_genome.rep3.fa LIU.rep3.tuxedo.bam.format.vcf.gz

#### pasa

singularity run docker://pasapipeline/pasapipeline
singularity shell pasapipeline_latest.sif

#### obsutil

obsutil share-cp ${TOKEN} ./  -ac=123456 -f -r  2>&1 | tee download.log

#### distmat

distmat -sequence total.aln -nucmethod 2 -outfile total.aln.distmat

#### busco

conda activate busco
cp *shortsummary* pasa.pep.busco.embryophyta.v4.1.2/
generate_plot.py -wd pasa.pep.busco.embryophyta.v4.1.2


#### OrthoFinder


orthofinder -b previous_orthofinder_directory -f new_fasta_directory

#### EDTA

CDS=input_maize/maize10_cds.fa
GENOME=input_maize/maize10.genome
# --species [Rice|Maize|others]
~/bin/EDTA/EDTA.pl  --species Maize --cds ${CDS} --curatedlib maizeTE02052020 --genome ${GENOME} --anno 1 --threads 40

#### parallel

ls *genome | parallel -j 5 samtools faidx {} 
cat list |parallel -j100 echo {}

#### fastp

#### MSScanX
i=abc.gff;python -m jcvi.formats.gff bed --type=mRNA --key=ID ${i} -o ${i%.gff3}.bed

MCScanX ./mcscan_run/AesCer

$ls mcscan_run/
tail -n +12 *ty |sed "s/^#.*/###/; s/.*:\s\+//" > ce_ae.anchors

#### gff_to_cds
mamba install gffread
gffread -x ${p}.cds -g ${p}.fa ${p}.gff   # ${p}.cds is output

#### conda environment
conda

##### yicluster
mamba install -c bioconda -c conda-forge  parallel exonerate phyx  bowtie2 bandage spades raxml-ng sra-tools phyx Beast   Gblocks corset Salmon TransDecoder  Trinity Bowtie2 FastQC Transrate  cd-hit Blat perl julia iqtree=2.1.2 MAFFT bioconductor-deseq2  r-ggplot2 bioconductor-ggtree tmux taxonkit ipython bioconductor-decipher r-mclust
mamba create -n eggnog -c bioconda noarch::eggnog-mapper=2.1.6
###### lht
mamba install gffread
mamba create -n r  r-ggplot2 bioconductor-ggtree ipython bioconductor-decipher bioconductor-deseq2 r-mclust
mamba create -n gatk4 python=3.6 gatk4

#### pip environment
pip install jcvi
iga
pip install gff3tool


#### Trinity

```
mamba create -n trinity -c bioconda cd-hit rsem trinity transdecoder

for i in `cat id2`
do
    bsub -o $i.log -e $i.err -n 30 -J trinity "Trinity --SS_lib_type RF         --seqType fq         --max_memory 50G --CPU 30         --left ${i}_1.clean.fq.gz         --right ${i}_2.clean.fq.gz         --output ${i}_trinity "
done

$TRINITY_HOME/util/TrinityStats.pl $i > $i.stat

$TRINITY_HOME/util/align_and_estimate_abundance.pl --seqType fq \
     --left data/GSNO_SRR1582648_1.fastq \
     --right data/GSNO_SRR1582648_2.fastq  \
     --transcripts trinity_out_dir/Trinity.fasta  \
     --output_prefix GSNO_SRR1582648 \
     --est_method RSEM  --aln_method bowtie2 \
     --trinity_mode --prep_reference \
     --output_dir GSNO_SRR1582648.RSEM
```

#### braker

```
source activate braker

export GENEMARK_PATH=/nfs/liuhui/bin/braker/dependency/gmes_linux_64
export PROTHINT_PATH=/nfs/liuhui/bin/braker/dependency/ProtHint/bin
export CDBTOOLS_PATH=/nfs/liuhui/bin/braker/dependency/cdbfasta
export PATH=/nfs/liuhui/bin/braker/dependency/GUSHR:$PATH
export MAKEHUB_PATH=/nfs/liuhui/bin/braker/dependency/MakeHub-1.0.5
export PATH=/nfs/liuhui/bin/braker/dependency:$PATH
export PATH=/nfs/liuhui/bin/braker/dependency/Augustus/bin:/nfs/liuhui/bin/braker/dependency/Augustus/scripts:$PATH
export AUGUSTUS_CONFIG_PATH=/nfs/liuhui/bin/braker/dependency/Augustus/config
export PATH=/nfs/liuhui/bin/braker/BRAKER/scripts:$PATH

module load boost

PEP=odb10_plants.fasta
THREADS=48

#for GENOME in Dipteryx_alata.fa.masked.fa
for GENOME in `cat sp_list`
do
    while [ ! -e ${GENOME} ]
    do
        echo "${GENOME} is not yet generated"
        sleep 1m
    done
    braker.pl -gff3 --cores=${THREADS} --genome=${GENOME} --prot_seq=$PEP --softmasking --workingdir braker_${GENOME}
done
```

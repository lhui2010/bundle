#!/bin/bash

# copied from Yuxing Xu

if [ -z $2 ]
then
    echo "Usage: gatk.sh fq1 fq2 reference(bwa-indexed) sample thread"
    exit 
fi


fq1=$1
fq2=$2
reference=$3
sample=$4
thread=$5

samtools=samtools
gatk=gatk

mkdir tmp
tmp_dir=`realpath tmp`
export TMPDIR=$tmp_dir

time bwa mem -t $thread -M -Y -R "@RG\tID:foo\tPL:ILLUMINA\tSM:$sample" \
    $reference $fq1 $fq2 | $samtools view -Sb - > ${sample}.bam && \
    echo "** BWA MEM done **"

time $samtools sort -@ $thread -m 4G -O bam -o ${sample}.sorted.bam ${sample}.bam && echo "** sorted raw bam file done **"

source activate gatk4

time $gatk MarkDuplicates -I ${sample}.sorted.bam -M ${sample}.markdup_metrics.txt -O ${sample}.sorted.markdup.bam && echo "** ${sample}.sorted.bam MarkDuplicates done **"

time $samtools index ${sample}.sorted.markdup.bam && echo "** ${sample}.sorted.markdup.bam index done **"

time $gatk HaplotypeCaller \
    -R $reference \
    --emit-ref-confidence GVCF \
    -I ${sample}.sorted.markdup.bam \
    -O ${sample}.g.vcf && echo "** gvcf done **"

time $gatk GenotypeGVCFs \
    -R $reference \
    -V ${sample}.g.vcf \
    -O ${sample}.vcf && echo "** vcf done **"

bgzip -f ${sample}.vcf
tabix -p vcf ${sample}.vcf.gz

time $gatk SelectVariants \
    -select-type SNP \
    -V ${sample}.vcf.gz \
    -O ${sample}.snp.vcf.gz

time $gatk VariantFiltration \
    -V ${sample}.snp.vcf.gz \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "PASS" \
    -O ${sample}.snp.filter.vcf.gz

time $gatk SelectVariants \
    -select-type INDEL \
    -V ${sample}.vcf.gz \
    -O ${sample}.indel.vcf.gz

time $gatk VariantFiltration \
    -V ${sample}.indel.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "PASS" \
    -O ${sample}.indel.filter.vcf.gz


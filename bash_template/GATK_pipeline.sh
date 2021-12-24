#!/bin/bash

#run_gatk.sh

#fastq1=../input/fastq1
#fastq2=../input/fastq2

REF="../input/bwa"

#step1 mapping
#bwa index $REF.fa
#bwa mem $REF.fa  -t 20  $fastq1 $fastq2 >bwa.sam 2>bwa.err


#step2 polishing

samtools view -bS $REF.sam >$REF.bam

java -jar ../bin/picard.jar CreateSequenceDictionary REFERENCE=$REF.fa OUTPUT=$REF.dict

java -Djava.io.tmpdir=`pwd`/tmp  -jar ../bin/picard.jar SortSam I=$REF.bam O=$REF.picard_sort.bam SO=coordinate TMP_DIR=`pwd`/tmp

java -Djava.io.tmpdir=`pwd`/tmp  -jar ../bin/picard.jar MarkDuplicates I=$REF.picard_sort.bam  O=$REF.picard_sort.dedup.bam METRICS_FILE=metrics.txt MAX_FILE_HANDLES=1000 TMP_DIR=`pwd`/tmp
java -Djava.io.tmpdir=`pwd`/tmp  -jar ../bin/picard.jar AddOrReplaceReadGroups I=$REF.picard_sort.dedup.bam  O=$REF.picard_sort.dedup.addrgp.bam  RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample1 TMP_DIR=`pwd`/tmp


BAM=$REF.picard_sort.dedup.addrgp.bam

#java -Xmx50g -jar ../bin/picard.jar BuildBamIndex INPUT=$BAM

#samtools faidx $REF

#if $$REF.fa=genome.fa faidx should be genome.fa.fai

#find relign targets
#java -jar ../bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF.fa -I $BAM  -o $BAM.realignment_targets.list

#realign
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R $REF.fa -I $BAM -targetIntervals $BAM.realignment_targets.list -o realigned_reads.bam

#call raw vcf
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF.fa -I realigned_reads.bam -o raw_variants.vcf

#extract snp & indel
java -jar GenomeAnalysisTK.jar -T SelectVariants -R $REF.fa -V raw_variants.vcf -selectType SNP -o raw_snps.vcf
java -jar GenomeAnalysisTK.jar -T SelectVariants -R $REF.fa -V raw_variants.vcf -selectType INDEL -o raw_indels.vcf


#filter snp
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R $REF.fa -V raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o filtered_snps.vcf


#filter indel
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R $REF.fa -V raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o filtered_indels.vcf


#Recalibration Score #1
java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF.fa -I realigned_reads.bam -knownSites filtered_snps.vcf -knownSites filtered_indels.vcf -o recal_data.table


#Recalibration Score #2
java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R $REF.fa -I realigned_reads.bam -knownSites filtered_snps.vcf -knownSites filtered_indels.vcf -BQSR recal_data.table -o post_recal_data.table


#summary of two recalibration
java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates -R $REF.fa -before recal_data.table -after post_recal_data.table -plots recalibration_plots.pdf


#Apply BQSR
java -jar GenomeAnalysisTK.jar -T PrintReads -R $REF.fa -I realigned_reads.bam -BQSR recal_data.table -o recal_reads.bam


#Call Variants #2
java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R $REF.fa -I recal_reads.bam -o raw_variants_recal.vcf


#Extract SNPs & Indels #2
java -jar GenomeAnalysisTK.jar -T SelectVariants -R $REF.fa -V raw_variants_recal.vcf -selectType SNP -o raw_snps_recal.vcf
java -jar GenomeAnalysisTK.jar -T SelectVariants -R $REF.fa -V raw_variants_recal.vcf -selectType INDEL -o raw_indels_recal.vcf


#Filter SNPs #2
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R $REF.fa -V raw_snps_recal.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o filtered_snps_final.vcf


#Filter Indels #2
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R $REF.fa -V raw_indels_recal.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o filtered_indels_recal.vcf



#java -Xmx50g -jar ../bin/GenomeAnalysisTK.jar -R $REF.fa -T HaplotypeCaller -I $BAM -o $REF.gatk.haplo.vcf
#java -Xmx50g -jar ../bin/GenomeAnalysisTK.jar -R $REF.fa -T FastaAlternateReferenceMaker -o $REF.gatk.haplo.fa --variant $REF.gatk.haplo.vcf

#mkdir ../output
#mv $REF.gatk.haplo.fa ../output/

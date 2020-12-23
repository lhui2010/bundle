

clean_mgiseq_sh="""
ADAPTER=/ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/adapters/MGISeq.fa
TAILCROP=145
HEADCROP=10
LEFT=${SAMPLE}_1.fq.gz
RIGHT=${SAMPLE}_2.fq.gz
java -jar  /ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/trimmomatic-0.38.jar PE \
-phred33 ${LEFT} ${RIGHT} \
${SAMPLE}_newadp_1.clean.fq.gz ${SAMPLE}_newadp_1.clean.unpair.fq.gz \
${SAMPLE}_newadp_2.clean.fq.gz ${SAMPLE}_newadp_2.clean.unpair.fq.gz \
ILLUMINACLIP:${ADAPTER}:0:30:10 \
LEADING:3 TRAILING:3 CROP:${TAILCROP} HEADCROP:${HEADCROP} SLIDINGWINDOW:1:10 MINLEN:75
"""
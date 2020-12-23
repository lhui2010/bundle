

clean_mgiseq_sh="""
ADAPTER=/ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/adapters/MGISeq.fa
TAILCROP=145
HEADCROP=10
LEFT={0}
RIGHT={1}
java -jar  /ds3200_1/users_root/yitingshuang/applications/Trimmomatic-0.38/trimmomatic-0.38.jar PE \
-phred33 $LEFT $RIGHT \
$LEFT.clean.fq.gz $LEFT.clean.unpair.fq.gz \
$RIGHT.clean.fq.gz $RIGHT.clean.unpair.fq.gz \
ILLUMINACLIP:${{ADAPTER}}:0:30:10 \
LEADING:3 TRAILING:3 CROP:$TAILCROP HEADCROP:$HEADCROP SLIDINGWINDOW:1:10 MINLEN:75
"""
"""
maker relevant utils
"""
import argparse
import os
import re
import sys
import time
from collections import defaultdict
import os.path as op

from parse import parse

from iga.apps.base import ActionDispatcher, sh, conda_act, workdir_sh, logger, Config, abspath_list, split_fasta, mkdir, \
    mv, wait_until_finish, bsub, emain

# def sam2gff(sam, gff=""):
#
# isoseq_official_sh=r"""
# ccs [movie].subreads.bam [movie].ccs.bam --min-rq 0.9
# lima --isoseq --dump-clips --no-pbi --peek-guess -j 24 ccs.bam primers.fasta demux.bam
# isoseq3 refine --require-polya combined_demux.consensusreadset.xml primers.fasta flnc.bam
# bamtools convert -format fastq -in flnc.bam > flnc.fastq
# """


# 0 subreads.bam
# 1 workdir
# 2 output.bam
isoseq_sh = r"""export PATH=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/isoseq3/BGI-Full-Length-RNA-Analysis-Pipeline/bin:$PATH
export PERL5LIB=""
WORKDIR={1}
mkdir ${{WORKDIR}}
ccs {0} ${{WORKDIR}}/ccs.bam --min-passes 0 --min-length 50 --max-length 21000 --min-rq 0.9
cd ${{WORKDIR}}
samtools view ccs.bam | awk '{{print ">"$1"\n"$10}}' > ccs.fa
echo ">primer_F
AAGCAGTGGTATCAACGCAGAGTACATGGGGGGGG
>primer_S
GTACTCTGCGTTGATACCACTGCTTACTAGT">primer.fa
makeblastdb -in primer.fa -dbtype nucl
blastn -num_threads 32 -query ccs.fa -db primer.fa -outfmt 7 -word_size 5 > mapped.m7
classify_by_primer -blastm7 mapped.m7 -ccsfa ccs.fa -umilen 8 -min_primerlen 16 -min_isolen 200 -outdir ./
flnc2sam ccs.sam isoseq_flnc.fasta > isoseq_flnc.sam
samtools view -bS isoseq_flnc.sam > isoseq_flnc.bam
isoseq3 cluster isoseq_flnc.bam unpolished.bam --verbose --use-qvs
#isoseq3 cluster isoseq_flnc.bam unpolished.bam --split-bam 10
# pbindex subreads.bam
# for i in `seq 0 9`
# do
#     isoseq3 polish unpolished.${{i}}.bam ${{ROOT}}/input/*.subreads.bam polished.${{i}}.bam --verbose &
# done
# wait
# samtools merge -@ 20 polished_total.bam polished.*.bam
# isoseq3 summarize polished_total.bam summary.csv
# Result is polished_total.bam.fastq 
"""


# TODO: allow multiple subreads input
def isoseq(subreads=None, workdir=''):
    r"""
    isoseq subreads.fasta

    Wrapper for `isoseq`
    """
    if (type(subreads) == list):
        subreads = " ".join(subreads)
    if (workdir == ''):
        workdir = "workdir_isoseq_" + subreads.split()[0]
    cmd = conda_act.format('isoseq3') + isoseq_sh.format(subreads, workdir)
    sh(cmd)


# 0 workdir
# 1 subreads.bam
# 2 primer.fa
isoseq_pb_sh = r"""mkdir -p {0}
cd {0}
ln -s ../{1}
ln -s ../{2}
ccs {1} {1}.ccs.bam --min-rq 0.9
#lima is used to remove primer sequence
#but can it be used to identify reads containing primer sequence as full length reads?
lima {1}.ccs.bam {2} {1}.fl.bam --isoseq  --peek-guess
#flnc equals full-length, non-concatemer
INPUTBAM={1}
PRIMER={2}
isoseq3 refine ${{INPUTBAM%.bam}}.fl.*.bam ${{PRIMER}} ${{INPUTBAM%.bam}}.flnc.bam --require-polya
isoseq3 cluster ${{INPUTBAM%.bam}}.flnc.bam ${{INPUTBAM%.bam}}.clustered.bam --verbose --use-qvs
"""


# TODO: allow multiple subreads input
def isoseq_pb(subreads=None, workdir=''):
    r"""
    convert Isoseq(pacbio standard) subreads.bam to flnc.fastq
    :param subreads:
    :param workdir:
    :return:
    """
    if (type(subreads) == list):
        subreads = " ".join(subreads)
    if (workdir == ''):
        workdir = "workdir_isoseq_" + subreads.split()[0]
    cmd = conda_act.format('isoseq3') + isoseq_sh.format(subreads, workdir)
    sh(cmd)


# 0 ref genome
# 1 qry fasta
fastq2gff_sh = """
minimap2 -t20 -C5 -ax splice {0} {1} |samtools view -F 256 -b - \
| bedtools bamtobed -split -i - > {1}.bed 
gt bed_to_gff3 {1}.bed | sort -k9,9 -k1,1 -k7,7 -k4,4n > {1}.rawgff
"""


def fastq2gff(fastq=None, genome=None, workdir=''):
    r"""
    align fastqs generated by Trinity or Isoseq to
    references and transform to maker acceptable gffs
    """
    # fastq = fastq)
    # genome = str(genome)
    logger.debug(fastq)
    logger.debug(genome)
    #    if (workdir == ''):
    #        workdir = "workdir_fastq2gff_" + fastq
    # workdir_sh.format(workdir) +
    cmd = conda_act.format('EDTA') + \
          fastq2gff_sh.format(genome, fastq)
    sh(cmd)
    rawgff = fastq + '.rawgff'
    gff = format_gt_gff_to_maker_gff(rawgff)
    return gff


def format_gt_gff_to_maker_gff(gff=None):
    """
    Format gff generated by gt to maker acceptable format
    Input Example:
    000034F|arrow_np1212    .       BED_feature     1899770 1900081 60      -       .       Name=TRINITY_GG_10000_c0_g1_i1
    000034F|arrow_np1212    .       BED_feature     2360009 2360227 60      -       .       Name=TRINITY_GG_10001_c0_g1_i1
    000034F|arrow_np1212    .       BED_feature     2359766 2359901 60      -       .       Name=TRINITY_GG_10001_c0_g2_i1
    000034F|arrow_np1212    .       BED_feature     2360103 2360189 60      -       .       Name=TRINITY_GG_10001_c0_g2_i1
    000034F|arrow_np1212    .       BED_feature     2360284 2360402 60      -       .       Name=TRINITY_GG_10001_c0_g2_i1
    Output Example:
    000034F|arrow_np1212    match   Trinity_Minimap 1899770 1900081 0       -       .       ID=LL_rep3-+-TRINITY_GG_10000_c0_g1_i1-+-0;Name=LL_rep3-+-TRINITY_GG_10000_c0_g1_i1-+-0
    000034F|arrow_np1212    match_part      Trinity_Minimap 1899770 1900081 60      -       .       ID=LL_rep3-+-TRINITY_GG_10000_c0_g1_i1-+-0-exon-1;Name=LL_rep3-+-TRINITY_GG_10000_c0_g1_i1-+-0;Parent=LL_rep3-+-TRINITY_GG_10000_c0_g1_i1-+-0
    000034F|arrow_np1212    match   Trinity_Minimap 2360009 2360227 0       -       .       ID=LL_rep3-+-TRINITY_GG_10001_c0_g1_i1-+-0;Name=LL_rep3-+-TRINITY_GG_10001_c0_g1_i1-+-0
    000034F|arrow_np1212    match_part      Trinity_Minimap 2360009 2360227 60      -       .       ID=LL_rep3-+-TRINITY_GG_10001_c0_g1_i1-+-0-exon-0;Name=LL_rep3-+-TRINITY_GG_10001_c0_g1_i1-+-0;Parent=LL_rep3-+-TRINITY_GG_10001_c0_g1_i1-+-0
    000034F|arrow_np1212    match   Trinity_Minimap 2359766 2361129 0       -       .       ID=LL_rep3-+-TRINITY_GG_10001_c0_g2_i1-+-0;Name=LL_rep3-+-TRINITY_GG_10001_c0_g2_i1-+-0
    """
    qry1_file = gff
    output_file = gff + ".gff"
    output_buff = ""
    prefix = re.sub(r'\..*', '', qry1_file)
    intron_cutoff = 20000
    source = 'Trinity_Minimap'
    last_chr = ''
    last_end = ''
    last_strand = ''
    last_pos = []
    last_lines = []
    last_name = ''
    count = 0
    feat_count = defaultdict(int)
    with open(qry1_file) as fh:
        for line in fh:
            if (line.startswith('#')):
                continue
            mylist = line.rstrip().split()
            # print(mylist[-1])
            mylist[-1] = mylist[-1].replace('/', '_')
            # Name feats
            this_chr = mylist[0]
            this_type = 'match_part'
            this_start = mylist[3]
            this_end = mylist[4]
            this_score = mylist[5]
            this_strand = mylist[6]
            this_phase = mylist[7]
            feat_parse = parse("Name={name}", mylist[-1])
            try:
                feat_name = prefix + "-+-" + feat_parse['name'] + "-+-" + str(
                    feat_count[feat_parse['name']])  # + this_chr.replace('|', '') + this_start
            except TypeError:
                logger.warning("TypeError on feat {}".format(mylist[-1]))
            this_feat = "ID={}-exon-{};Name={};Parent={}".format(feat_name, str(count), feat_name, feat_name)
            line = "\t".join(
                [this_chr, this_type, source, this_start, this_end, this_score, this_strand, this_phase, this_feat])

            if (last_chr == '' or
                    last_name == feat_name and
                    last_chr == this_chr and
                    last_strand == this_strand and
                    abs(last_end - int(this_start)) < intron_cutoff):
                count += 1

                this_feat = "ID={}-exon-{};Name={};Parent={}".format(feat_name, str(count), feat_name, feat_name)
                line = "\t".join(
                    [this_chr, this_type, source, this_start, this_end, this_score, this_strand, this_phase, this_feat])
                last_chr = this_chr
                last_end = int(this_end)
                last_pos.append(int(this_start))
                last_pos.append(int(this_end))
                last_strand = this_strand
                last_lines.append(line)
                last_name = feat_name

            else:
                if (last_name == feat_name and last_chr == this_chr and last_strand == this_strand):
                    logger.warning("LargeIntron {} on {}".format(str(abs(last_end - int(this_start))), feat_name))
                # Prepare print
                match_chr = last_chr
                match_type = "match"
                match_start = str(min(last_pos))
                match_end = str(max(last_pos))
                match_score = '0'
                match_strand = last_strand
                match_phase = '.'
                match_feat = "ID={};Name={}".format(last_name, last_name)
                match_line = "\t".join(
                    [match_chr, match_type, source, match_start, match_end, match_score, match_strand, match_phase,
                     match_feat])
                last_lines.insert(0, match_line)
                # print("\n".join(last_lines))
                output_buff += "\n".join(last_lines) + "\n"
                # Initialize
                count = 0
                last_chr = this_chr
                last_end = int(this_end)
                last_pos = [int(this_start), int(this_end)]
                last_strand = this_strand
                # Add count of same name so no duplicate name would occur
                feat_name_slim = last_name.split('-+-')[1]
                # exit()
                feat_count[feat_name_slim] += 1
                # print(feat_name_slim)
                # print(feat_count[feat_name_slim])
                try:
                    feat_name = prefix + "-+-" + feat_parse['name'] + "-+-" + str(feat_count[feat_parse['name']])
                except TypeError:
                    logger.warning("TypeError on feat {}".format(mylist[-1]))
                # print(feat_parse['name'])
                # print(feat_count[feat_parse['name']])
                this_feat = "ID={}-exon-{};Name={};Parent={}".format(feat_name, str(count), feat_name, feat_name)
                line = "\t".join(
                    [this_chr, this_type, source, this_start, this_end, this_score, this_strand, this_phase, this_feat])
                last_lines = [line]
                last_name = feat_name
        else:
            match_chr = last_chr
            match_type = "match"
            match_start = str(min(last_pos))
            match_end = str(max(last_pos))
            match_score = '0'
            match_strand = last_strand
            match_phase = '.'
            match_feat = "ID={};Name={}".format(last_name, last_name)
            match_line = "\t".join(
                [match_chr, match_type, source, match_start, match_end, match_score, match_strand, match_phase,
                 match_feat])
            last_lines.insert(0, match_line)
            # print("\n".join(last_lines))
            output_buff += "\n".join(last_lines) + "\n"
    with open(output_file, 'w') as fh:
        fh.write(output_buff)
    return output_file


# environment for maker
maker_env_sh = r"""
ZOE_HMM_DIR=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/snap/Zoe/HMM/

"""

# 0 workdir
# 1 prev_round
# 2 current_round
# 3 cfg_file
# 4 cfg

maker_run_sh = r"""
cd {}
maker *ctl >> maker.out 2>> maker.err
"""


def maker_run(genome=None, estgff=None, pepgff=None,
                 rmgff=None, round=1, species='', use_grid='T',
                 augustus_species='', snap_hmm=''):
    """
    Give genome and evidence, run maker gene prediction in parallel
    :param genome:
    :param estgff:
    :param pepgff:
    :param rmgff:
    :param round:
    :param species:
    :param use_grid: whether to use LSF to submit jobs
    :param augustus_species: species for augustus, if
        from round1, then it's directory name, like coriaria_contig.fa_R1
    :param snap_hmm: hmm file for snap, if from round1,
        then it's directory name, like coriaria_contig.fa_R1
    :return:
    """
    snap_hmm_dir = '/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/snap/Zoe/HMM/'
    workdir = ''
    if (species == ''):
        workdir = genome + '_R' + str(round)
    #logger.warning(workdir)
    #exit(1)
    # Split genome and placing into working directory like:
    # coriaria_round1:
    #   chunk.1/1.fa
    #   chunk.2/2.fa
    fa_list = split_fasta(genome, workdir, 100)
    # default returned a string with file names, changing it into list type
    # logger.debug(fa_list)
    # change estgff file name in to absolute path
    # logger.warning(os.getcwd())
    [estgff, pepgff, rmgff] = abspath_list([estgff, pepgff, rmgff])
    #    logger.warning([estgff, pepgff, rmgff])
    # Preparing cfg files
    cfg_exe = Config('maker_exe')
    cfg_bopts = Config('maker_bopts')
    cfg = Config('maker')
    cfg.update('est_gff={};protein_gff={};rm_gff={}'.format(estgff, pepgff, rmgff))
    if (round == 1):
        #Only first round will be ran in direct predict mode
        cfg.update('est2genome=1;protein2genome=1')
    else:
        cfg.update('est2genome=0;protein2genome=0')
        if('.hmm' not in snap_hmm):
            #In case .hmm extension was not added in input
            snap_hmm = snap_hmm + '.hmm'
        cfg.update('snaphmm={0};augustus_species={1}'.format(
            op.join(snap_hmm_dir, snap_hmm), augustus_species))
    #
    # get abs path of all fasta files
    os.chdir(workdir)
    fa_list = abspath_list(fa_list)
    # abspath_list(fa_list)
    # job list for storing submitted job IDs
    job_list = []
    for i in fa_list:
        fa_name = op.basename(i)
        workdir_sep = i + '.run/'
        mkdir(workdir_sep)
        mv(i, workdir_sep)
        # fasta = op.join(workdir, fa_name)
        cfg.update('genome={}'.format(fa_name))
        cfg.write_to_file(op.join(workdir_sep, "maker_opts.ctl"))
        cfg_exe.write_to_file(op.join(workdir_sep, 'maker_exe.ctl'))
        cfg_bopts.write_to_file(op.join(workdir_sep, 'maker_bopts.ctl'))
        cmd = maker_run_sh.format(workdir_sep, cfg)
        # sh(cmd)
        job_id = bsub(cmd)
        job_list.append(job_id)

    logger.warning("Submitted jobs:")
    logger.warning(job_list)
    wait_until_finish(job_list)
    logger.warning("Submmited job finished, check log files to make sure they really finished")


def maker_check(workdir=None):
    """
    check whether all partitions finished
    :param workdir:
    :return:
    """
    subdir = os.listdir(workdir)
    error_list = []
    unfinished_list = []

    for sd in subdir:
        if ('run' in sd):
            maker_log = op.join(workdir, sd, 'maker.err')
            maker_log_buff = ''
            fail_mark = 1
            try:
                with open(maker_log) as fh:
                    maker_log_buff = fh.read()
                if ('Maker is now finished!!!' in maker_log_buff):
                    if (not 'ERROR' in maker_log_buff and not 'Fail' in maker_log_buff):
                        pass
                    else:
                        error_list.append(sd)
                else:
                    unfinished_list.append(sd)
            except FileNotFoundError:
                logger.error("Can't find maker error log for {}".format(sd))
                unfinished_list.append(sd)
    if (len(unfinished_list) + len(error_list) == 0):
        logger.warning("Cheers! All finished without errors!")
    else:
        if (len(unfinished_list) > 0):
            logger.warning("Unfinished chunks are:")
            [logger.warning(l) for l in unfinished_list]
        if (len(error_list) > 0):
            logger.warning("Chunks with errors are:")
            [logger.warning(l) for l in error_list]
        exit(1)

    return error_list + unfinished_list


def deploy_augustus():
    r"""
    deploy augustus config dir to /tmp/lh/ for maker use
    :return:
    """
    node_list = {}
    node_list['Q64C1T_X4'] = ['node02', 'node03', 'node04', 'node05']
    node_list['Q104C512G_X4'] = ['node10', 'node11', 'node12', 'node13']
    augustus_config_dir = '/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/config'
    local_dir = '/tmp/lh_config'
    cmd = 'touch {0} && rm -fr {0} && cp -fr {1} {0}'.format(local_dir, augustus_config_dir)
    for q in node_list:
        for node in node_list[q]:
            bsub(cmd, q + ' -m {}'.format(node))
    time.sleep(120)
    return 0

# 0 working directory
collect_maker_sh = r"""
cd {}
#9.fa.run/9.maker.output/9_master_datastore_index.log
#000041F|arrow_np1212    9_datastore/6C/BE/000041F%7Carrow_np1212/       STARTED
#000041F|arrow_np1212    9_datastore/6C/BE/000041F%7Carrow_np1212/       FINISHED
touch total_master_datastore_index.log
for i in `ls -d *.fa.run/`
do
    echo $i
    j=${{i%.fa.run/}}
    cat ${{i}}/${{j}}.maker.output/${{j}}_master_datastore_index.log |sed "s/\t/\t${{j}}.fa.run\/${{j}}.maker.output\//" \
    >>total_master_datastore_index.log
done


#b. Merge fasta
fasta_merge -d total_master_datastore_index.log
gff3_merge -o genome.all.gff -d total_master_datastore_index.log
gff3_merge -n -o genome.all.noseq.gff -d total_master_datastore_index.log

echo "Merge completed succefully:"

awk '$2=="maker"' genome.all.noseq.gff > genome.maker.gff
# echo "##FASTA" >> genome.maker.gff
# cat *.run/*.fa >> genome.maker.gff
cat *.run/*.fa > ref.fa

date"""


def maker_collect(workdir=None):
    """
    Collect maker result from a paralleled run in workdir
    :param workdir:
    :return:
    """
    cmd = collect_maker_sh.format(workdir)
    res = sh(cmd)
    logger.warning(res)
    return 0


busco_export_sh = r"""
export BUSCO_CONFIG_FILE=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/busco/myconfig.ini
"""

# 0 workdir
# 1 PREFIX of this model
train_snap_sh = r"""export HMMDIR=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/snap/Zoe/HMM/
mkdir -p {0}/train_snap
cd {0}/train_snap
if [ ! -e genome.all.gff ]
then
    ln -s ../genome.all.gff
fi
maker2zff -x 0.25 -l 50  genome.all.gff
fathom -gene-stats genome.ann genome.dna >gene-stats.log 2>&1
fathom -validate genome.ann genome.dna >validate.log 2>&1

fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna

mkdir -p params

cd params

forge ../export.ann ../export.dna >../forge.log 2>&1

cd ..
hmm-assembler.pl snap_trained params > snap_trained.hmm

#if [ -d ${{HMMDIR}}/{1}.hmm ]
#then
#    RND=$(date +%s%N)
#    mv ${{HMMDIR}}/{1}.hmm ${{HMMDIR}}/{1}.hmm.$RND
#fi

cp snap_trained.hmm ${{HMMDIR}}/{1}.hmm

echo "Train SNAP completed succefully:"
echo "${{HMMDIR}}/{1}.hmm"
date
"""

# 0 workdir
# 1 prefix
train_augustus_sh = r"""
mkdir -p {0}/train_augustus
cd {0}/train_augustus
if [ ! -e genome.all.gff ]
then
    ln -s ../genome.all.gff
fi

if [ ! -e ref.fa ]
then
    ln -s ../ref.fa
fi

awk -v OFS="\t" '{{ if ($3 == "mRNA") print $1, $4, $5 }}' genome.all.gff | \
  awk -v OFS="\t" '{{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }}' | \
  bedtools getfasta -fi ref.fa -bed - -fo total.all.maker.transcripts1000.fasta

cp  /ds3200_1/users_root/yitingshuang/lh/projects/buzzo/maker/../busco/myconfig.ini  ./config.ini
export BUSCO_CONFIG_FILE=$PWD/config.ini
#If not export to local dir, PREFIX_exon_probs.pbl will not be generated and maker will fail
export AUGUSTUS_CONFIG_PATH=/tmp/lh_config

LINEAGE=viridiplantae_odb10
LINEAGE=embryophyta_odb10
THREADS=104
INPUT=total.all.maker.transcripts1000.fasta
OUTPUT={1}
NEWMODEL={1}
#TODO
AUGUSTUS_SPECIES=arabidopsis
AUGUSTUS_CONFIG_PATH_ORIGINAL=/ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/config
if [ -d $OUTPUT ]
then
    rm -rf $OUTPUT
fi
busco -i $INPUT  -o $OUTPUT  -l $LINEAGE \
  -m genome -c $THREADS --long --augustus_species $AUGUSTUS_SPECIES \
  --augustus_parameters='--progress=true' >busco.out 2>busco.err

cd $OUTPUT/run_${{LINEAGE}}/augustus_output/retraining_parameters/BUSCO_${{OUTPUT}}/

rename "BUSCO_" "" *

sed -i 's/BUSCO_//g' {1}_parameters.cfg

if [ -d $AUGUSTUS_CONFIG_PATH_ORIGINAL/species/$NEWMODEL ]
then
    RND=$(date +%s%N)
    mv $AUGUSTUS_CONFIG_PATH_ORIGINAL/species/$NEWMODEL $AUGUSTUS_CONFIG_PATH_ORIGINAL/species/$NEWMODEL.$RND
    #mv $AUGUSTUS_CONFIG_PATH/species/$NEWMODEL $AUGUSTUS_CONFIG_PATH/species/${{NEWMODEL}}.$RND
fi
mkdir -p $AUGUSTUS_CONFIG_PATH_ORIGINAL/species/$NEWMODEL
cp ./${{OUTPUT}}*  $AUGUSTUS_CONFIG_PATH_ORIGINAL/species/{1}/
#
echo "Train Augustus completed succefully"
echo "Augustus species: {1}"
date

"""


def maker_train(workdir=None, prefix='', augustus='T', snap='T', use_grid='T'):
    """
    :param workdir:
    :return:
    """
    cmd = ''
    workdir = op.abspath(workdir)
    set_workdir = 'cd {};'.format(workdir)
    if (prefix == ''):
        prefix = op.basename(workdir)
    if (snap == 'T'):
        cmd += set_workdir + "\n" + train_snap_sh.format(workdir, prefix)
    if (augustus == 'T'):
        #BUSCO 4.1.2 failed to retrain augustus, even specified augustus config dir to local
        #The error was no exon_probs.pbl file produced.
        #I don't know why. For now, I used busco v4.0.1(with bug manual fixed).
        #cmd += set_workdir + "\n" + conda_act.format('busco') + train_augustus_sh.format(workdir, prefix)
        cmd += set_workdir + "\n" + busco_export_sh + train_augustus_sh.format(workdir, prefix)
    if (use_grid == 'T'):
        joblist = bsub(cmd, direct_submit='F')
        wait_until_finish(joblist)
    else:
        sh(cmd)
    return 0


def str_to_class(str1):
    return getattr(sys.modules[__name__], str1)


# def minimap_rna(transcript, genome, threads=30, output=''):
#     if (output == ''):
#         output = transcript + ".sam"
#     cmd = minimap_rna_sh.format(threads, genome, transcript, output)
#     sh(cmd)


# parallel run

# train

# liftover
# Require RaGOO


def liftover():
    pass


# function

def main():
    """
    the main function
    """
    emain()
    # actions = (
    #     ('isoseq', 'extract isoseq flnc reads from subreads.bam')
    #     ('fastq2gff', 'map fastq to reference genome and get gff files'),
    # )


    # p = ActionDispatcher(actions)
    # p.dispatch(globals())

    # print(__file__)
    # print(__doc__)
    # exit()
    # prog_name = "busco_wrapper"
    # usage = "run busco on selected GENOME"
    #
    # parser = argparse.ArgumentParser(
    #     prog=prog_name,
    #     formatter_class=argparse.RawDescriptionHelpFormatter,
    #     description=textwrap.dedent(usage),
    #     epilog="")
    # parser.add_argument("GENOME", help="Genome to be evalutated in fasta format")
    # parser.add_argument("-t", "--threads", default=64, type=int, help="flanking distance default (1000)")
    # args = parser.parse_args()
    #
    # busco(args.GENOME)


#    flanking_distance = args.flanking

if __name__ == "__main__":
    main()

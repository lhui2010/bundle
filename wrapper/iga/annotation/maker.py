"""
maker relevant utils
"""
import os
import os.path as op
import re
import sys
import time
from collections import defaultdict

from parse import parse

from iga.apps.base import sh, conda_act, logger, Config, abspath_list, split_fasta, mkdir, \
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


def fix_comma_in_parent():
    r"""
    change Parent=CORNET00004845-1,CORNET00004845-2 to
    ...Parent=CORNET00004845-1
    ...Parent=CORNET00004845-2
    Read from stdin and print to stdout
    :return:
    """
    import fileinput
    fileinput.sys.argv = sys.argv[2:]
    for l in fileinput.input():
        found_error = re.search(r'(.*;)Parent=(.*,.*)', l)
        if (found_error):
            leading = found_error[1]
            lagging = found_error[2].rstrip(';').split(',')
            for lag in lagging:
                print(leading + "Parent=" + lag + ';')
        else:
            print(l, end='')


# 0 reference_genome.fa
# 1 transcripts.fa
# 2 maker.gff
pasa_refine_sh = r"""

export PASAHOME=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/maker/bin/PASApipeline.v2.4.1
export PATH=$PASAHOME/bin/:$PATH
export PERL5LIB=$PASAHOME/PerlLib:$PASAHOME/SAMPLE_HOOKS:/ds3200_1/users_root/yitingshuang/lh/anaconda3/lib/site_perl/5.26.2/

touch {0}.sqlite

echo "## templated variables to be replaced exist as <__var_name__>

# database settings
DATABASE=$PWD/{0}.sqlite

#######################################################
# Parameters to specify to specific scripts in pipeline
# create a key = "script_name" + ":" + "parameter"
# assign a value as done above.

#script validate_alignments_in_db.dbi
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=80
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=80

#script subcluster_builder.dbi
subcluster_builder.dbi:-m=50
" >   pasa.alignAssembly.sqlite.txt

#-+- Align transcript to genome and create the original sqlite db
$PASAHOME/Launch_PASA_pipeline.pl \
   -c pasa.alignAssembly.sqlite.txt -C -R -g {0} \
   -t {1}  \
    --ALIGNERS gmap,blat --CPU 20
  
#-+- Rename maker.gff(necessary?)

#-+- Sort maker.gff in gene, mRNA, exon, CDS order (other type is not needed by pasa)
mv {2} {2}.bak
awk '$3=="gene"' {2}.bak > {2}
awk '$3=="mRNA"' {2}.bak >> {2}
awk '$3=="CDS"' {2}.bak >> {2}
#-+- Fix the Parent=mRNA1,mRNA2 issue. changing them into two lines. Maker is so weird :(
awk '$3=="exon"' {2}.bak |python -m iga.annotation.maker fix_comma_in_parent >> {2}

#-+- Load GFF
$PASAHOME/scripts/Load_Current_Gene_Annotations.dbi \
    -c pasa.alignAssembly.sqlite.txt \
    -g {0} \
    -P {2}
    
echo "# database settings
DATABASE=$PWD/{0}.sqlite

#######################################################
# Parameters to specify to specific scripts in pipeline
# create a key = script_name + : + parameter
# assign a value as done above.


#script cDNA_annotation_comparer.dbi
cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP=<__MIN_PERCENT_OVERLAP__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_PROT_CODING=<__MIN_PERCENT_PROT_CODING__>
cDNA_annotation_comparer.dbi:--MIN_PERID_PROT_COMPARE=<__MIN_PERID_PROT_COMPARE__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_FL_COMPARE=<__MIN_PERCENT_LENGTH_FL_COMPARE__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_NONFL_COMPARE=<__MIN_PERCENT_LENGTH_NONFL_COMPARE__>
cDNA_annotation_comparer.dbi:--MIN_FL_ORF_SIZE=<__MIN_FL_ORF_SIZE__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_ALIGN_LENGTH=<__MIN_PERCENT_ALIGN_LENGTH__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP_GENE_REPLACE=<__MIN_PERCENT_OVERLAP_GENE_REPLACE__>
cDNA_annotation_comparer.dbi:--STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE=<__STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE__>
cDNA_annotation_comparer.dbi:--TRUST_FL_STATUS=<__TRUST_FL_STATUS__>
cDNA_annotation_comparer.dbi:--MAX_UTR_EXONS=<__MAX_UTR_EXONS__>
cDNA_annotation_comparer.dbi:--GENETIC_CODE=<__GENETIC_CODE__>
" > pasa.annotCompare.config

#-+- Refine genome gff, (Very slow)
$PASAHOME/Launch_PASA_pipeline.pl \
    -c pasa.annotCompare.config -A \
    -g {0} \
    -t {1}

"""


def pasa_refine(genome=None, transcript=None, gff=None, use_grid='F'):
    r"""
    :param genome: the assembled genome (fasta)
    :param transcript: the assembled transcripts (fasta)
    :param gff: the existing annotation (gff3)
    :param use_grid: whether to use bsub (T or F)
    :return:
    """
    cmd = pasa_refine_sh.format(genome, transcript, gff)
    if(use_grid == 'T'):
        bsub(cmd, direct_submit='F', cpus=5)
    else:
        sh(cmd)


# 0 input.gff
# 1 prefix, like A188 or CORNE
maker_rename_sh = r"""
maker_map_ids --abrv_gene 'G' --prefix {1} --justify 8 --suffix '-t' --iterate 1 {0}  > {0}.map.txt
#Caution map_gff_ids will rewrite the file instead of generating a new one
cp {0} {0}.format.gff
map_gff_ids {0}.map.txt {0}.format.gff
"""

def maker_rename_gff(gff=None, prefix='MAKER'):
    cmd = maker_rename_sh.format(gff, prefix)
    sh(cmd)

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
              rmgff=None, round=1, species='', use_grid='T', cpus=1,
              augustus_species='', snap_hmm='', queue='Q104C512G_X4', update=''):
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
    if (os.path.exists(workdir)):
        rnd = str(time.time())
        mv(workdir, workdir + rnd)
    # logger.warning(workdir)
    # exit(1)
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
        # Only first round will be ran in direct predict mode
        cfg.update('est2genome=1;protein2genome=1')
    else:
        cfg.update('est2genome=0;protein2genome=0')
        if ('.hmm' not in snap_hmm):
            # In case .hmm extension was not added in input
            snap_hmm = snap_hmm + '.hmm'
        cfg.update('snaphmm={0};augustus_species={1}'.format(
            op.join(snap_hmm_dir, snap_hmm), augustus_species))
    if (update != ''):
        cfg.update(update)
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
        workdir_sep = op.abspath(workdir_sep)
        # fasta = op.join(workdir, fa_name)
        cfg.update('genome={}'.format(fa_name))
        cfg.write_to_file(op.join(workdir_sep, "maker_opts.ctl"))
        cfg_exe.write_to_file(op.join(workdir_sep, 'maker_exe.ctl'))
        cfg_bopts.write_to_file(op.join(workdir_sep, 'maker_bopts.ctl'))
        cmd = maker_run_sh.format(workdir_sep)
        # sh(cmd)
        if (use_grid == 'T'):
            job_id = bsub(cmd, queue=queue, cpus=2)
            job_list.append(job_id)
            time.sleep(30)
        else:
            job_list.append(cmd)
    if (use_grid == 'T'):
        logger.warning("Submitted jobs:")
        logger.warning(job_list)
        wait_until_finish(job_list)
        logger.warning("Submmited job finished, check log files to make sure they really finished")
    else:
        sh(job_list, parallel='T', cpus=cpus)


maker_resub_sh = r"""
cd {}
mkdir -p rm 
mv *.maker.output rm/
rm -rf rm &
maker *ctl > maker.out 2> maker.err
"""


def maker_resub(dir_list=None, queue="Q104C512G_X4", cpus=4):
    r"""
    Resubmit failed jobs by directory name
    :param dir_list:
    :param queue:
    :return:
    """
    if (type(dir_list) == str):
        dir_list = [dir_list]
    logger.warning(dir_list)
    # logger.debug(queue)
    # exit(1)
    job_list = []
    for i in dir_list:
        cmd = maker_resub_sh.format(i)
        job_id = bsub(cmd, queue=queue, cpus=cpus)
        job_list.append(job_id)
        time.sleep(30)
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
        # exit(1)

    return error_list + unfinished_list


def maker_check_resub(workdir=None, queue="Q104C512G_X4"):
    i = 1
    current_dir = op.abspath(os.curdir)
    # absworkdir = op.abspath(workdir)
    while (i < 3):
        i += 1
        # at most resub two times
        os.chdir(current_dir)
        failed_list = maker_check(workdir)
        os.chdir(workdir)
        if (len(failed_list) > 1):
            maker_resub(failed_list, queue=queue, cpus=i)
        else:
            return 0
    logger.error("Failed too many times, stop submitting")
    exit(1)
    return 1


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
set -euxo pipefail
cd {}
#9.fa.run/9.maker.output/9_master_datastore_index.log
#000041F|arrow_np1212    9_datastore/6C/BE/000041F%7Carrow_np1212/       STARTED
#000041F|arrow_np1212    9_datastore/6C/BE/000041F%7Carrow_np1212/       FINISHED
touch total_master_datastore_index.log
rm    total_master_datastore_index.log
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

# added optimise step as BUSCO will fail this step naturally, maybe need 5 hours to finish
# Error: IO.c: loadable library and perl binaries are mismatched (got handshake key 0xdb80080, needed 0xde00080)
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

#Do not need it in current HPC environment
#export AUGUSTUS_CONFIG_PATH=/tmp/lh_config

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
fi
mkdir -p $AUGUSTUS_CONFIG_PATH_ORIGINAL/species/$NEWMODEL
cp ./${{OUTPUT}}*  $AUGUSTUS_CONFIG_PATH_ORIGINAL/species/{1}/

#/ds3200_1/users_root/yitingshuang/lh/anaconda3/envs/busco/scripts/optimize_augustus.pl --cpus=8 \
#--species={1} ../../training_set.db
#
echo "Train Augustus completed succefully"
echo "Augustus species: {1}"
date

"""


def filter_gff_by_aed(gff=None, gff_out='', aed='0.2'):
    r"""
    filter maker_gff by aed value
    output as gff.filter as default
    :param gff:
    :param aed:
    :return:
    """
    buff_list = []
    buff = ''
    if (gff_out == ''):
        gff_out = gff + '.filter'
    with open(gff) as fh:
        for line in fh:
            buff += line
            mylist = line.split()
            if (mylist[2] == 'gene'):
                buff_list.append(buff)
                buff = ''
    result = ''
    for bf in buff_list:
        pattern_result = re.search(r'_AED=(.*?);', bf)[1]
        if (float(pattern_result) <= float(aed)):
            result += bf
    with open(gff_out) as fh:
        fh.write(result)
    return 0


# training augustus without BUSCO
# run after snap is finished
# single/8 thread, default in local run.
# 0 workdir
# 1 prefix
# 2 absolute path to full length fasta
train_augustus_direct_sh = r"""
if [ -d {0}/train_augustus_direct ]
then
    rm -rf {0}/train_augustus_direct
fi
mkdir -p {0}/train_augustus_direct
cd {0}/train_augustus_direct
if [ ! -e genome.all.gff ]
then
    ln -s ../genome.all.gff
fi

if [ ! -e ref.fa ]
then
    ln -s ../ref.fa
fi

NUMFOUND=500
NUMSPLIT=250
CDNA_FASTA={2}
AUGUSTUS_SPECIES_NAME={1}_direct
WORKING_DIR=$PWD
ROOT=$PWD

if [ -d $AUGUSTUS_CONFIG_PATH/species/$AUGUSTUS_SPECIES_NAME ]
then
    RND=$(date +%s%N)
    mv $AUGUSTUS_CONFIG_PATH/species/$AUGUSTUS_SPECIES_NAME $AUGUSTUS_CONFIG_PATH/species/$AUGUSTUS_SPECIES_NAME.$RND
fi

ln -s ../train_snap/uni.ann
ln -s ../train_snap/uni.dna 

/ds3200_1/users_root/yitingshuang/lh/bin/GC_specific_MAKER/fathom_to_genbank.pl --annotation_file uni.ann \
 --dna_file uni.dna  --genbank_file augustus.gb \
 --number ${{NUMFOUND}}
perl -e  'while (my $line = <>){{ if ($line =~ /^LOCUS\s+(\S+)/) {{ print "$1\n"; }} }}'  ${{WORKING_DIR}}/augustus.gb \
 >  ${{WORKING_DIR}}/genbank_gene_list.txt
/ds3200_1/users_root/yitingshuang/lh/bin/GC_specific_MAKER/get_subset_of_fastas.pl  \
 -l  ${{WORKING_DIR}}/genbank_gene_list.txt   \
 -f ${{WORKING_DIR}}/uni.dna  -o  ${{WORKING_DIR}}/genbank_gene_seqs.fasta
 
perl ~/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/scripts/randomSplit.pl ${{WORKING_DIR}}/augustus.gb ${{NUMSPLIT}}

~/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/scripts/autoAug.pl --species=$AUGUSTUS_SPECIES_NAME \
--genome=${{WORKING_DIR}}/genbank_gene_seqs.fasta --trainingset=${{WORKING_DIR}}/augustus.gb --cdna=$CDNA_FASTA  \
--noutr

# Failed due to mysql issue
# --pasa --useGMAPforPASA

cd ./autoAug/autoAugPred_abinitio/shells

x=1
while [ -e ./aug${{x}} ]
do
    echo "A.  $x"
    ./aug${{x}} &
    let x=x+1
done

wait

cd $WORKING_DIR

~/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/scripts/autoAug.pl --species=$AUGUSTUS_SPECIES_NAME \
--genome=${{WORKING_DIR}}/genbank_gene_seqs.fasta --useexisting --hints=${{WORKING_DIR}}/autoAug/hints/hints.E.gff \
 -v -v -v  --index=1

cd ${{WORKING_DIR}}/autoAug/autoAugPred_hints/shells/

let x=1
while [ -e ./aug${{x}} ]
do
    echo "B.  $x"
    ./aug${{x}} &
    let x=x+1
done

wait

echo "Successfully finished"
"""


def maker_train(workdir=None, prefix='', augustus='T', snap='T', use_grid='T', augustus_direct='T',
                cdna_fasta=''):
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
        # BUSCO 4.1.2 failed to retrain augustus, even specified augustus config dir to local
        # The error was no exon_probs.pbl file produced.
        # I don't know why. For now, I used busco v4.0.1(with bug manual fixed).
        # cmd += set_workdir + "\n" + busco_export_sh + train_augustus_sh.format(workdir, prefix)
        cmd += set_workdir + "\n" + conda_act.format('busco') + busco_export_sh + \
               train_augustus_sh.format(workdir, prefix)
    if (augustus_direct == 'T'):
        if (cdna_fasta != ''):
            cdna_fasta = op.abspath(cdna_fasta)
            logger.warning(workdir)
            logger.warning(cdna_fasta)
            #TODO in future: use pasa to train augustus, may need to add pasa_export_sh before this line
            cmd += train_augustus_direct_sh.format(workdir, prefix, cdna_fasta)
        else:
            logger.error("Provide cdna.fasta before train augustus_direct")
            exit(1)
    if (use_grid == 'T'):
        joblist = bsub(cmd, direct_submit='F', cpus=2)
        wait_until_finish(joblist)
    else:
        sh(cmd)
    return 0


def str_to_class(str1):
    return getattr(sys.modules[__name__], str1)



# liftover
# Require RaGOO
liftover_sh = r"""
lift_over.py 
"""

def liftover_by_agp(gff=None, agp=None):
    """
    lift over gff files based on agp file, this is used when annotation genes on contig level and need to transfer to
    chromosome level
    :param gff:
    :param agp:
    :return:
    eg:
    python -m iga.annotation.maker liftover_by_agp genome.maker.gff CORNE_v1.0.agp > genome.lifted.gff
    gff_genome_to_genes.pl genome.lifted.gff.cds CORNE_v1.0.chr.fa > genome.lifted.gff.cds.fa
    cds2aa.pl genome.lifted.gff.cds.fa > genome.lifted.gff.cds.fa.pep
    """
    #chr name of contig
    reverse_strand = {'-':'+', '+':'-'}
    chrd = {}
    #start on chr
    startd = {}
    #end on chr
    endd = {}
    #contig strand on chr
    strandd = {}
    #contig length
    lengthd = {}
    with open(agp) as fh:
        for line in fh:
            if (line.startswith('#')):
                continue
            else:
                mylist = line.split()
                contig_name = mylist[5]
                if(mylist[4] == 'W'):
                    chrd[contig_name]    = mylist[0]
                    startd[contig_name]  = int(mylist[1])
                    endd[contig_name]    = int(mylist[2])
                    strandd[contig_name] = mylist[-1]
                    lengthd[contig_name] = int(mylist[-2])

    with open(gff) as fh:
        for line in fh:
            if(line.startswith('#')):
                print(line, end='')
            else:
                mylist = line.rstrip().split()
                this_contig = mylist[0]
                this_start = int(mylist[3])
                this_end = int(mylist[4])
                this_strand = mylist[6]
                #transforming
                new_chr = chrd[this_contig]
                if(strandd[this_contig] == '-'):
                    new_strand = reverse_strand[this_strand]
                    #new start coordinate = seqLength - endCoord
                    this_start = lengthd[this_contig] - (this_start - 1)
                    this_end = lengthd[this_contig] - (this_end - 1)
                else:
                    new_strand = this_strand
                #Add offset, because loci in AGP is 1-based, so minus 1 is the real offset
                new_start = this_start + startd[this_contig] - 1
                new_end = this_end + startd[this_contig] - 1
                #in case we need swap loci like Chr0 . gene 800 1
                if(new_end < new_start):
                    (new_start, new_end) = (new_end, new_start)
                mylist[0] = new_chr
                mylist[3] = str(new_start)
                mylist[4] = str(new_end)
                mylist[6] = new_strand
                new_line = "\t".join(mylist)
                print(new_line)



# function

# def main():
#     """
#     the main function
#     """
#     emain()
#     # actions = (
#     #     ('isoseq', 'extract isoseq flnc reads from subreads.bam')
#     #     ('fastq2gff', 'map fastq to reference genome and get gff files'),
#     # )
#
#     # p = ActionDispatcher(actions)
#     # p.dispatch(globals())
#
#     # print(__file__)
#     # print(__doc__)
#     # exit()
#     # prog_name = "busco_wrapper"
#     # usage = "run busco on selected GENOME"
#     #
#     # parser = argparse.ArgumentParser(
#     #     prog=prog_name,
#     #     formatter_class=argparse.RawDescriptionHelpFormatter,
#     #     description=textwrap.dedent(usage),
#     #     epilog="")
#     # parser.add_argument("GENOME", help="Genome to be evalutated in fasta format")
#     # parser.add_argument("-t", "--threads", default=64, type=int, help="flanking distance default (1000)")
#     # args = parser.parse_args()
#     #
#     # busco(args.GENOME)


#    flanking_distance = args.flanking

if __name__ == "__main__":
    emain()

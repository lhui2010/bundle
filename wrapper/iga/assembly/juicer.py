'''
hic relevant scripts
'''


from iga.apps.base import bsub, get_prefix, abspath_list
import os.path as op


# 0 prefix
# 1 genome.fa
# 2 hic.fastq
# 3 threads
# 4 enzyme
juicer_pipe_sh = '''
ROOT=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/juicer
PREFIX={0}
WORKDIR=workdir_juicer_{0}
JUICERDIR=$ROOT/juicer
THREADS={3}
GENOME={1}

touch ${{WORKDIR}}
rm -r ${{WORKDIR}}
mkdir -p ${{WORKDIR}}


cd ${{WORKDIR}}
ln -s ${{ROOT}}/juicer/misc
ln -s ${{ROOT}}/juicer/scripts

mkdir -p ${{WORKDIR}}/fastq
cd ${{WORKDIR}}/fastq

for i in {2}
do 
    ln -s $i
done

mkdir -p ${{WORKDIR}}/references
cd ${{WORKDIR}}/references

ln -s ${{GENOME}}

REF=`basename ${{GENOME}}`

cd ${{WORKDIR}}

bwa index ${{REF}}
samtools faidx ${{REF}}
cut -f1,2 ${{REF}}.fai > ${{REF}}.len
FALEN=${{REF}}.len

mkdir -p  ${{WORKDIR}}/restriction_sites

python misc/generate_site_positions.py {4} ${{PREFIX}} ${{GENOME}}

#${PREFIX}_MboI.txt
ENZYME_FILE=${{PREFIX}}_{4}.txt
mv ${{ENZYME_FILE}} ${{WORKDIR}}/restriction_sites/
cd ${{WORKDIR}}
${{ROOT}}/juicer/scripts/juicer.sh.mnd_only  -z ${{REF}} -p ${{FALEN}} \\
  -y ${{WORKDIR}}/restriction_sites/${{ENZYME_FILE}} -D ${{WORKDIR}} -d ${{WORKDIR}} -t ${{THREADS}} \\
   >${{WORKDIR}}/jc.out 2>${{WORKDIR}}/jc.err

mkdir -p ${{WORKDIR}}/3ddna
cd ${{WORKDIR}}/3ddna
bash ${{ROOT}}/3d-dna/run-asm-pipeline.sh  -m diploid ${{GENOME}} ${{WORKDIR}}/aligned/merged_nodups.txt

'''


def juicer_pipe(genome=None, hic_fastq=None, prefix='', cpus=40, enzyme='MboI', queue='Q104C512G_X4'):
    """
    :param genome: genome.fasta
    :param hic_fastq:  "hic1.fastq hic2.fastq"
    :param prefix: prefix of your species
    :param enzyme: [MboI|]
    :param cpus: threads, default 40
    :param queue: Q104C512G_X4
    :return:
    """
    # 0 prefix
    # 1 genome.fa
    # 2 hic.fastq
    # 3 threads
    # 4 enzyme
    if prefix == '':
        prefix = get_prefix(genome)
    genome = op.abspath(genome)
    fq_lists = hic_fastq.split()
    abspath_list(fq_lists)
    hic_fastq = '"'
    for fq in fq_lists:
        hic_fastq += fq + " "
    hic_fastq = hic_fastq.rstrip() + '"'
    cmd = juicer_pipe_sh.format(prefix, genome, hic_fastq, cpus, enzyme)
    bsub(cmd, cpus=cpus, name='juicer.'+prefix, queue=queue)
#Input:
#1. Genome file
#2. HiC fastq file

#Output:
#1. .hic file

def get_hic():
    """ Use juicer and 3d-dna to obtain mnd and 0.hic file
    """

def post():
    """ Use .mnd and curated .assembly file to get .agp and chr.fasta file
    """

def assembly2agp():
    """ converts assembly file to agp file
    """

def liftover():
    """ lift over gff annotation from contig to chromosome 
    """

def buildchr():
    """ build chromosome based on .assembly or .agp
    """

'''
hic relevant scripts
'''

# 0 prefix
# 1 genome.fa
# 2 hic.fastq
# 3 threads
'''
ROOT=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/juicer
PREFIX={0}
WORKDIR=workdir_juicer_{0}
touch ${{WORKDIR}}
rm -r ${{WORKDIR}}
mkdir -p ${{WORKDIR}}
# ln -s ../input_fastq ${{WORKDIR}}/fastq
JUICERDIR=$ROOT/juicer
THREADS={3}
GENOME={1}

bwa index ${{GENOME}}
samtools faidx ${{GENOME}}
cut -f1,2 ${{GENOME}}.fai > ${{GENOME}}.len
FALEN=${{GENOME}}.len
#mkdir fastq
#mv <fastq_files> fastq/
#juicer/scripts/juicer.sh -D  ${JUICERDIR} -d ${WORKDIR} -t ${THREADS}

mkdir -p ${{WORKDIR}}/references ${{WORKDIR}}/restriction_sites
cd ${{WORKDIR}}
ln -s ${{ROOT}}/juicer/misc
ln -s ${{ROOT}}/juicer/scripts
python misc/generate_site_positions.py MboI ${{PREFIX}} ${{GENOME}}

#${PREFIX}_MboI.txt
ENZYME_FILE=${{PREFIX}}_MboI.txt
mv ${ENZYME_FILE} ${WORKDIR}/restriction_sites/
cd ${WORKDIR}
${{ROOT}}/juicer/scripts/juicer.sh.mnd_only  -z ${{GENOME}} -p ${{FALEN}} \\
  -y ${{WORKDIR}}/restriction_sites/${{ENZYME_FILE}} -D ${{WORKDIR}} -d ${{WORKDIR}} -t ${{THREADS}} \\
   >${{WORKDIR}}/jc.out 2>${WORKDIR}/jc.err


mkdir -p ${WORKDIR}/3ddna
cd ${WORKDIR}/3ddna
bash ${ROOT}/3d-dna/run-asm-pipeline.sh  -m diploid ${GENOME} ${WORKDIR}/aligned/merged_nodups.txt
'''

def run_juicer():
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

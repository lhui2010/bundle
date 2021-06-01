'''
hic relevant scripts
'''
import logging
import re
from collections import defaultdict, OrderedDict

from iga.apps.base import bsub, get_prefix, abspath_list, emain
import os.path as op


class AssemblyIO:
    def __init__(self, assembly_file):
        """
        Input assembly file
        #==> allhic.0.review.assembly <==
        #>Hic.fastq.gz.counts_GATC.20g1:::fragment_1 1 5423057
        #>Hic.fastq.gz.counts_GATC.20g1:::fragment_2:::debris 2 50000
        #>Hic.fastq.gz.counts_GATC.20g1:::fragment_3 3 7785000
        #...
        #1 -3 13 22 -5
        """
        chr_id = 0
        # Storing the contig order information that be used in jcbat grouping, eg dic[1]=contig1
        self.order_to_frag_name = {}
        # Storing the contig size information that be used in jcbat grouping, eg dic[1]=contig1
        self.order_to_frag_size = {}
        # chr_dict storing the grouping of
        self.chr_dict = defaultdict(list)
        self.chr_size = defaultdict(int)
        with open(assembly_file) as fh:
            for line in fh:
                (this_chr, this_start, this_end, this_name, this_strand) = ["."] * 5
                # Each line is a new chromosome
                if line.startswith(">"):
                    (frag_name, frag_order, frag_len) = line.rstrip().split()
                    frag_name = frag_name[1:]
                    frag_len = int(frag_len)
                    self.order_to_frag_name[frag_order] = frag_name
                    self.order_to_frag_size[frag_order] = frag_len
                else:
                    this_chr = "chr" + "%02d" % chr_id
                    mylist = line.rstrip().split()
                    for i in mylist:
                        self.chr_dict[this_chr].append(i)
                        self.chr_size[this_chr] += self.order_to_frag_size[i.replace('-', '')]
                    chr_id += 1

    def sort_by_size(self):
        """
        sort reviewed assembly file into large to small order
        :return:
        """
        chr_count = len(self.chr_size.keys())
        new_chr_dict = OrderedDict()
        new_chr_size = OrderedDict()
        final_chr_list = []
        sorted_chr = sorted(self.chr_size, key=self.chr_size.__getitem__, reverse=True)
        for i in range(0, chr_count):
            chr_id = i + 1
            this_chr = "chr" + "%02d" % chr_id
            final_chr_list.append(this_chr)
            new_chr_dict[this_chr] = self.chr_dict[sorted_chr[i]]
            new_chr_size[this_chr] = self.chr_size[sorted_chr[i]]
        self.chr_dict = new_chr_dict
        self.chr_size = new_chr_size

    def gettext(self, size=False):
        """
        :return: return the text format of current assembly file
        """
        output = ''
        for k in self.order_to_frag_name:
            output += ">{} {} {}\n".format(self.order_to_frag_name[k], k, self.order_to_frag_size[k])
        for k in self.chr_dict:
            if size:
                output += str(self.chr_size[k]) + ":"
            output += " ".join(self.chr_dict[k])
            output += "\n"
        return output

    def toagp(self):
        """
        :return: return the text of current assembly file in AGP format
        """
        pass


# 0 prefix
# 1 genome.fa
# 2 hic.fastq
# 3 threads
# 4 enzyme
juicer_pipe_sh = r'''
ROOT=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/juicer
PREFIX={0}
WORKDIR=$PWD/workdir_juicer_{0}
JUICERDIR=$ROOT/juicer
THREADS={3}
GENOME={1}

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

mv *1.*gz HiC_R1.fastq.gz
mv *2.*gz HiC_R2.fastq.gz

mkdir -p ${{WORKDIR}}/references
cd ${{WORKDIR}}/references

ln -s ${{GENOME}}

REF=`basename ${{GENOME}}`

bwa index ${{REF}}
samtools faidx ${{REF}}
cut -f1,2 ${{REF}}.fai > ${{REF}}.len
FALEN=${{REF}}.len

cd ${{WORKDIR}}
mkdir -p  ${{WORKDIR}}/restriction_sites
python misc/generate_site_positions.py {4} ${{PREFIX}} ${{GENOME}}

ENZYME_FILE=${{PREFIX}}_{4}.txt
mv ${{ENZYME_FILE}} ${{WORKDIR}}/restriction_sites/
cd ${{WORKDIR}}
${{ROOT}}/juicer/scripts/juicer.sh.mnd_only  -z ${{WORKDIR}}/references/${{REF}} \
  -p ${{WORKDIR}}/references/${{FALEN}} \
  -y ${{WORKDIR}}/restriction_sites/${{ENZYME_FILE}} -D ${{WORKDIR}} -d ${{WORKDIR}} -t ${{THREADS}} \
   >${{WORKDIR}}/jc.out 2>${{WORKDIR}}/jc.err

mkdir -p ${{WORKDIR}}/3ddna
cd ${{WORKDIR}}/3ddna
bash ${{ROOT}}/3d-dna/run-asm-pipeline.sh  -m diploid ${{GENOME}} ${{WORKDIR}}/aligned/merged_nodups.txt

'''


def juicer_pipe(genome=None, hic_fastq=None, prefix='', threads=40, enzyme='MboI', queue='Q104C512G_X4',
                submit='T'):
    """
    :param genome: genome.fasta
    :param hic_fastq:  "hic1.fastq hic2.fastq"
    :param prefix: prefix of your species
    :param enzyme: [MboI|HindIII]
    :param threads: threads, default 40
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
    hic_fastq = ''
    for fq in fq_lists:
        hic_fastq += fq + " "
    hic_fastq = hic_fastq.rstrip()
    cmd = juicer_pipe_sh.format(prefix, genome, hic_fastq, threads, enzyme)
    if submit == 'T':
        bsub(cmd, cpus=threads, name='juicer.' + prefix, queue=queue)


# Input:
# 1. Genome file
# 2. HiC fastq file

# Output:
# 1. .hic file


def get_hic():
    """ Use juicer and 3d-dna to obtain mnd and 0.hic file
    """
    return 0


post_assembly_sh = r"""
ROOT=$PWD
BIN=/ds3200_1/users_root/yitingshuang/lh/projects/buzzo/juicer/3d-dna/run-asm-pipeline-post-review.sh

ASSEMBLY={0}
CONTIG={1}
MND={2}

${{BIN}} -r ${{ASSEMBLY}} ${{CONTIG}} ${{MND}}
"""


def post_assembly(assembly=None, contig=None, mnd_file=None):
    """ Use .mnd and curated .assembly file to get .agp and chr.fasta file
    :param assembly: The review.assembly file generated by juicerbox
    :param contig: The contig fasta file
    :param mnd_file: The merged non-duplicated file, usually can be found in aligned/merged_nodups.txt of juicer dir
    """
    cmd = post_assembly_sh.format(assembly, contig, mnd_file)
    bsub(cmd, name="post_assembly")
    logging.warning("""Example Output:
falcon_v340_sgs_polish.FINAL.assembly: contig order information, should be the same with reviewed.assembly
falcon_v340_sgs_polish.FINAL.fasta: Assembled pseudochromosomes
falcon_v340_sgs_polish.final.fasta: Contig file in fasta format
falcon_v340_sgs_polish.final.hic:  .HiC file to be reviewed in juicerbox
""")
    return 0


def sort_assembly(assembly=None, size='F'):
    """
    sort reviewed assembly in large to small order
    :param assembly:
    :return:
    """
    Asb = AssemblyIO(assembly)
    Asb.sort_by_size()
    if size == 'T':
        print(Asb.gettext(size=True))
    else:
        print(Asb.gettext())


def assembly2agp():
    """ converts assembly file to agp file
    """
    ##assembly
    # >tig00000064|arrow_np1212 1 7029267
    # >tig00004986|arrow_np1212 2 6940838
    # >tig00004984|arrow_np1212 3 6851042
    # >tig00005003|arrow_np1212 4 6793917
    # ...
    # 235 -144 245 -150 -193 602 -1209 1211 505 131 157 -154 -85 127 280 80
    # 160
    # 1095
    # 1349
    # 859

    ##AGP
    # Hic.fastq.gz.counts_GATC.20g10 1   132119  1   W   000083F|arrow_np1212    1   132119  -
    # Hic.fastq.gz.counts_GATC.20g10 132120  132219  2   U   100 contig  yes map
    # Hic.fastq.gz.counts_GATC.20g10 132220  266408  3   W   000093F|arrow_np1212    1   134189  +
    # Hic.fastq.gz.counts_GATC.20g10 266409  266508  4   U   100 contig  yes map

    return 0


def liftover():
    """ lift over gff annotation from contig to chromosome 
    """
    return 0


def buildchr():
    """ build chromosome based on .assembly or .agp
    """
    return 0


if __name__ == '__main__':
    emain()

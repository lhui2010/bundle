#!/usr/bin/env python

import argparse
import textwrap
import pysam


#def sam2clipped_fq():

def get_qual_str(input_list):
    new_str = ''
    for l in input_list:
        new_str += chr(l + 33)
    return new_str

    


def main():
    prog_name = "Another python program"
    usage = "Another python program"

    parser = argparse.ArgumentParser(
        prog=prog_name, 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("qry1", help="qry1 file")
#    parser.add_argument("-f", "--flanking", default=10000, type=int, help="flanking distance default (1000)")
    args = parser.parse_args()  

    qry1_file = args.qry1
    samfile = pysam.AlignmentFile(qry1_file, "r")
#    flanking_distance = args.flanking
#    with open(qry1_file) as fh:
#        for line in fh:
#            mylist = line.rstrip().split()
    count=1
    len_cutoff = 75
    #'cigar', 'cigarstring', 'cigartuples'
    """
    +-----+--------------+-----+
    |M    |BAM_CMATCH    |0    |
    +-----+--------------+-----+
    |I    |BAM_CINS      |1    |
    +-----+--------------+-----+
    |D    |BAM_CDEL      |2    |
    +-----+--------------+-----+
    |N    |BAM_CREF_SKIP |3    |
    +-----+--------------+-----+
    |S    |BAM_CSOFT_CLIP|4    |
    +-----+--------------+-----+
    |H    |BAM_CHARD_CLIP|5    |
    +-----+--------------+-----+
    |P    |BAM_CPAD      |6    |
    +-----+--------------+-----+
    |=    |BAM_CEQUAL    |7    |
    +-----+--------------+-----+
    |X    |BAM_CDIFF     |8    |
    +-----+--------------+-----+
    |B    |BAM_CBACK     |9    |
    +-----+--------------+-----+
    |NM   |NM tag        |10   |
    +-----+--------------+-----+
    """
    
    with open("hic_clip_R1.fastq", 'w') as file1, open("hic_clip_R2.fastq", 'w') as file2:
        file_list = [file1, file2]
        read_tag_list = [r"\1", r"\2"]
#       print(read.cigar) [(4, 67), (0, 83)]
#       print(read.cigarstring) 67S83M
#       print(read.cigartuples) [(4, 67), (0, 83)]
#       print(read.get_cigar_stats()) (array('I', [83, 0, 0, 0, 67, 0, 0, 0, 0, 0, 0]), array('I', [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]))
#k        start = 0
#k        end = 150
#k        if('S' in read.cigarstring):
#k            if(read.cigar[0][0] == 4):
#k                start = read.cigar[0][1] - 1
#k            elif(read.cigar[-1][0] == 4):
#k                end = end - read.cigar[-1][1]
#               V300046408L2C001R0010000002
#                array('B', [32, 33, 34, 14, 36, 34, 32, 32, 37, 6
#                7, 37, 36, 36, 37, 34, 37, 37, 37, 19, 37, 37, 30
#                0, 30, 27, 37, 37, 37, 37, 37, 31, 37, 26, 37, 36
#                CTCTCTCTCTCTTACATGAACTTCAATCTGATACACAGTGAAGATGTTG
#                1179
#                CTCTCTCTCTCTTACATGAACTTCAATCTGATACACAGTGAAGATGTTG
#                array('B', [34, 21, 37, 37, 37, 37, 32, 37, 35, 3
#                37, 37, 16, 35, 37, 37, 37, 38, 37, 37, 38, 37, 3
#                34, 35, 35, 37, 36, 30, 38, 14, 34, 37, 37, 32, 1
#                GATCTGTGACGGCCGTAGGAGAATAGCTTGGAAGACTTTGTCTTTATGC
        for read in samfile.fetch():
            if(read.query_alignment_length <= len_cutoff):
                continue
            if(read.is_read1):
                index = 0
            if(read.is_read2):
                index = 1
            if(read.cigarstring):
#cigar string is not none
                file_list[index].write('@{}{}\n{}\n+\n{}\n'.format(read.query_name, 
                read_tag_list[index],
                read.query_alignment_sequence, 
                get_qual_str(read.query_alignment_qualities)))
            else:
#No S, then print entire region
                file_list[index].write('@{}{}\n{}\n+\n{}\n'.format(read.query_name, 
                read_tag_list[index],
                read.query_sequence, 
                get_qual_str(read.query_qualities)))
#            count +=1
#            if (count>10):
#                break

#        print(read.query_name)
#        print(read.cigarstring)
#        print(read.query_sequence)
#        print(read.query_qualities)
#        print(read.query_alignment_sequence)
#        print(read.query_alignment_qualities)
#        print(get_qual_str(read.query_alignment_qualities))
#        print(read.rname)
#        print(read.seq)
#        print(read.get_forward_qualities())
#        print(read.get_forward_sequence())
#        print(dir(read))
#        print(help(read))
#        if read.is_paired:
#            pairedreads.write(read)


if __name__ == "__main__":
    main()

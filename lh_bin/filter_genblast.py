#!/usr/bin/env python

import argparse
import textwrap
import re
import logging


class Feat:
    def __init__(self, gene_id, content, rank=-1, score=0, length=0, exon_num=0):
        self.gene_id = gene_id
        self.content = content
        self.exon_num = exon_num
        self.length = length
        rank_search = re.search(r'(.*)-R(\d+)', gene_id)
        if rank_search != None:
            self.core_gene_id = rank_search[1]
            self.rank = int(rank_search[2])
        if score == '0.0' or score == '0' or '-' in score or score == '.':
            self.score = 1
        else:
            self.score = float(score)

    def append(self, line):
        self.content += line
        self.exon_num += 1


class Best:
    def __init__(self, score, rank, gene_id):
        self.score = score
        self.rank = rank
        self.gene_id = gene_id


class GFF:
    def __init__(self, filename, transcript_key="transcript"):
        self.GFF_dict = {}
        self.best = {}
        with open(filename) as fh:
            transcript_id = ""
            score = 0
            count = 0
            for line in fh:
                count += 1
                if line.startswith('#') or line.rstrip() == "":
                    continue
                mylist = line.rstrip().split()
                if len(mylist) < 6 or len(mylist) < 3:
                    logging.warning("Error: Array element number is {} on \
                            line {}".format(len(mylist), count))
                    exit()
                if mylist[2] == transcript_key:
                    transcript_id = self.get_ID_from_last_feat(mylist[-1])
                    score = mylist[5]
                    gene_len = int(mylist[4]) - int(mylist[3])
                    self.GFF_dict[transcript_id] = \
                        Feat(gene_id=transcript_id, score=mylist[5], content=line, length=gene_len)

                    core_gene_id = self.GFF_dict[transcript_id].core_gene_id
                    this_score = self.GFF_dict[transcript_id].score
                    this_rank = self.GFF_dict[transcript_id].rank

                    if (core_gene_id not in self.best or
                            this_score > self.best[core_gene_id].score):
                        self.best[core_gene_id] = Best(this_score, this_rank, transcript_id)
                else:
                    #                    print(transcript_id)
                    #                    print(line)
                    self.GFF_dict[transcript_id].append(line)

    @staticmethod
    def get_ID_from_last_feat(lastfeat):
        lastfeat = re.sub(r";.*", "", lastfeat.replace('ID=', ''))
        return lastfeat

    def filter(self, div_threshold=0.9, rank_keep=2, exon_num=2, len_cutoff=50000):
        my_new_dict = {}
        for gene_id in self.GFF_dict:
            core_gene_id = self.GFF_dict[gene_id].core_gene_id
            try:
                div_result = self.GFF_dict[gene_id].score / self.best[core_gene_id].score
            except ZeroDivisionError as e:
                print("\t".join([gene_id, str(self.GFF_dict[gene_id].score), str(self.best[core_gene_id].score)]))

            if (self.GFF_dict[gene_id].rank > rank_keep or
                    self.GFF_dict[gene_id].score / self.best[core_gene_id].score <
                    div_threshold or
                    self.GFF_dict[gene_id].length > len_cutoff and
                    self.GFF_dict[gene_id].exon_num == exon_num or
                    self.GFF_dict[gene_id].length > len_cutoff and
                    self.GFF_dict[gene_id].rank > rank_keep):
                pass
            else:
                my_new_dict[gene_id] = self.GFF_dict[gene_id]

        self.GFF_dict = my_new_dict.copy()

    def print_out(self):
        for k in self.GFF_dict:
            print(self.GFF_dict[k].content, end='')


def main():
    prog_name = "Another python program"
    usage = "Another python program"

    parser = argparse.ArgumentParser(
        prog=prog_name,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage),
        epilog="")
    parser.add_argument("qry1", help="qry1 file")
    #    parser.add_argument("qry2", help="qry2 file")
    #    parser.add_argument("-f", "--flanking", default=10000, type=int, help="flanking distance default (1000)")
    args = parser.parse_args()

    qry1_file = args.qry1
    #    qry2_file = args.qry2
    #    flanking_distance = args.flanking

    #   000028F|arrow_np1212    genBlastG       transcript      2151627 2152061 23.3019 -       .       ID=sp_A7M942_PSAC_CUSGR
    #   000028F|arrow_np1212    genBlastG       coding_exon     2151627 2152061 .       -       .       ID=sp_A7M942_PSAC_CUSGR
    #   000005F|arrow_np1212    genBlastG       transcript      3483885 3484160 23.2295 -       .       ID=sp_A7M942_PSAC_CUSGR
    gff_read = GFF(qry1_file)
    gff_read.filter()
    gff_read.print_out()


if __name__ == "__main__":
    main()

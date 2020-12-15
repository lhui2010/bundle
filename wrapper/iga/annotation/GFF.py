"""
GFF relevant utils
"""
from collections import OrderedDict, defaultdict
from itertools import chain

from parse import parse

from iga.apps.base import emain, logger


class Feat:
    r"""
    The feat data structure that is needed by GFF class, support:
    1. Get parent
    2. Get childs
    #chr01   .       gene    12132486        12138762        .       -       .       ID=CORNEG00007591;
    #Name=CORNE00007591-t5;Alias=maker-000023F|arrow_np1212-snap-gene-26.41
    """

    def __init__(self, gff_line):
        self.childs = []
        self.content = gff_line.rstrip()
        mylist = self.content.split('\t')
        # Assign gff values by
        # https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        self.seqid = mylist[0]
        self.source = mylist[1]
        self.type = mylist[2]
        self.start = mylist[3]
        self.end = mylist[4]
        self.score = mylist[5]
        self.strand = mylist[6]
        self.phase = mylist[7]
        # Refactor affected vairables
        self.len = abs(int(self.end) - int(self.start)) + 1
        self.attributes = mylist[8]
        self.attr_dict = OrderedDict()
        attr_list = self.attributes.rstrip(';').split(';')
        for a in attr_list:
            try:
                (attr_key, attr_value) = parse("{}={}", a)
            except TypeError:
                logger.error("Type Error on line {}, list {} and attribute: {}".format(self.attributes, attr_list, a))
                continue
            self.attr_dict[attr_key] = attr_value
        if 'Parent' in self.attr_dict:
            self.parent = self.attr_dict['Parent']
        else:
            self.parent = None
        try:
            self.ID = self.attr_dict['ID']
        except KeyError:
            self.ID = None

    def get_parent(self):
        r"""
        Return the name of the parent
        :return:
        """
        return self.parent

    def add_child(self, child):
        r"""
        child is feat type
        :param child:
        :return:
        """
        if type(child) != Feat:
            raise TypeError
        self.childs.append(child)

    def get_all_child_feats(self, type=''):
        r"""
        return all lines that are descendants of this feat
        :param type:
        TODO: There i a bug here, will print gene mRNA CDS if use get_all_child_feats(CDS)
        :return:
        """
        if len(self.childs) == 0:
            if type == '' or self.type == type:
                # return when type is wild card or self.type equals specified type
                result = self.content + "\n"
        else:
            result = self.content + "\n"
            for i in self.childs:
                result += i.get_all_child_feats(type)
        return result

    def get_all_child_feats_obj(self, type=''):
        r"""
        return all descendants feats of this feat
        :param type:
        :return:
        """
        if type == '' or self.type == type:
            # return when type is wild card or self.type equals specified type
            result = [self]
        else:
            result = []
        if len(self.childs) == 0:
            pass
        else:
            for i in self.childs:
                result += i.get_all_child_feats_obj(type)
#        if any(isinstance(i, list) for i in result):
#            result = list(chain.from_iterable(result))
        return result

    def update_tag(self, tag, value):
        r"""
        Add or append a tag to this feat
        :param tag:
        :param value:
        :return:
        """
        self.attr_dict[tag] = value
        self.__refactor__()

    def delete_tag(self, tag):
        r"""
        Remove a tag
        :param tag:
        :return:
        """
        if tag == "ID":
            logger.error("Can't delete ID item")
            return 1
        del self.attr_dict[tag]
        self.__refactor__()
        return 0

    def __refactor__(self):
        r"""
        Update all fields
        :return:
        """
        self.attributes = ''
        for i in self.attr_dict:
            self.attributes += "{}={};".format(i, self.attr_dict[i])
        self.content = "\t".join([self.seqid, self.source, self.type, self.start, self.end,
                                  self.score, self.strand, self.phase, self.attributes])
        self.ID = self.attr_dict['ID']
        if 'Parent' in self.attr_dict:
            self.parent = self.attr_dict['Parent']
        return 0

    def print_all_childs(self):
        result = self.get_all_child_feats()
        print(result, end="")


class GFF:
    r"""
    GFF class that support:
    1. reading a GFF into memory
    2. adding feat by gene ID or transcript ID
    3. print out
    ...3. extracting feat by level or by transcript ID into tab delimited file
    """

    def __init__(self, filename):
        # Top level, which has no parent, usually gene type
        self.top_level_list = []
        # A dictionary store all type
        self.GFF_dict = OrderedDict()
        # cds name is same
        count_cds = defaultdict(int)
        with open(filename) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                feat = Feat(line)
                if "Parent" not in feat.content:
                    # Top level
                    self.top_level_list.append(feat.attr_dict['ID'])
                    self.GFF_dict[feat.ID] = feat
                else:
                    if feat.type == "CDS":
                        # Manage duplicate CDS
                        original_cds_name = feat.attr_dict['ID']
                        count_cds[original_cds_name] += 1
                        if count_cds[original_cds_name] > 1:
                            new_name = "{}:cds{}".format(original_cds_name, count_cds[original_cds_name])
                            feat.update_tag("ID", new_name)
                            if original_cds_name in self.GFF_dict:
                                new_first_cds = "{}:cds{}".format(original_cds_name, 1)
                                self.GFF_dict[new_first_cds] = self.GFF_dict[original_cds_name]
                                self.GFF_dict[new_first_cds].update_tag("ID", new_first_cds)
                                del self.GFF_dict[original_cds_name]
                            self.GFF_dict[new_name] = feat
                        else:
                            self.GFF_dict[feat.ID] = feat
                    else:
                        self.GFF_dict[feat.ID] = feat
                    parent = feat.parent
                    self.GFF_dict[parent].add_child(feat)

    def print_out(self):
        r"""
        print out gff to screen
        :return:
        """
        total_result = ''
        for k in self.top_level_list:
            result = self.GFF_dict[k].get_all_child_feats()
            print(result.rstrip())
        return 0

    def to_str(self):
        r"""
        print out gff to screen
        :return:
        """
        total_result = ''
        for k in self.top_level_list:
            result = self.GFF_dict[k].get_all_child_feats()
            total_result += result.rstrip() + "\n"
        return total_result

    def get_attr(self, attr=''):
        r"""
        Return attr as a dict, like extracting _AED from maker GFF
        :return:
        """
        result = OrderedDict()
        for k in self.GFF_dict:
            if attr in self.GFF_dict[k].attr_dict:
                result[k] = self.GFF_dict[k].attr_dict[attr]
        return result

    def longest_mRNA(self):
        longest_table = ""
        longest_gff = ""
        for k in self.top_level_list:
            mRNA_list = self.GFF_dict[k].get_all_child_feats_obj('mRNA')
            longest = 0
            if mRNA_list == []:
                # This is not a protein coding gene
                continue
            else:
                for mRNA in mRNA_list:
                    mRNA_len = 0
                    CDS_list = mRNA.get_all_child_feats_obj('CDS')
                    if not CDS_list:
                        logger.error("ERROR: mRNA {} do not have CDS type".format(mRNA))
                    else:
                        for CDS in CDS_list:
                            mRNA_len += CDS.len
                    mRNA.abs_len = mRNA_len
                    if longest < mRNA.abs_len:
                        longest = mRNA.abs_len
                        self.GFF_dict[k].longest = mRNA.ID
            longest_table += "{}\t{}".format(k, self.GFF_dict[k].longest) + "\n"
            longest_gff += self.GFF_dict[k].content + "\n"
            longest_gff += self.GFF_dict[self.GFF_dict[k].longest].get_all_child_feats()
        return [longest_table, longest_gff]


def longest_mRNA(gff=None):
    r"""
    print mRNA's longest to gff.longest.table and gff.longest.gff
    :param gff:
    :return:
    """
    gff_obj = GFF(gff)
    (longest_table, longest_gff) = gff_obj.longest_mRNA()
    with open(gff + "longest.table", 'w') as fh:
        fh.write(longest_table)
    with open(gff + "longest.gff", 'w') as fh:
        fh.write(longest_gff)
    return 0


def extract_gff_tag(gff=None, tag=None):
    r"""
    Extract specific attribute of gene or mRNA or any Item by the attribute name
    like python %s maker.gff _AED
    :param GFF:
    :param tag:
    :return:
    """
    gff_db = GFF(gff)
    tag_dict = gff_db.get_attr(tag)
    for k in tag_dict:
        print("{}\t{}".format(k, tag_dict[k]))

def fix_gt_gff(gff=None):
    r"""
    Fix the resulting gff files from gt gtf2gff3
    :param gff:
    :return:
    """
    count = defaultdict(int)
    with open(gff) as fh:
        for line in fh:
            if(line.startswith('#')):
                print(line, end ='')
                continue
            feat = Feat(line)
            if feat.type == "gene":
                feat.update_tag("ID",
                                feat.attr_dict['gene_id'].replace('gene:', ''))
            elif feat.type == "mRNA":
                feat.update_tag("Parent",
                                feat.attr_dict['gene_id'].replace('gene:', ''))
                feat.update_tag("ID",
                                feat.attr_dict['transcript_id'].replace('mRNA:', ''))
            else:
                feat.parent = feat.attr_dict['transcript_id'].replace('mRNA:', '')
                prefix = "{}:{}".format(feat.parent + feat.type)
                count[prefix] += 1
                feat.update_tag("ID",
                                "{}-{}".format(prefix, count[prefix]))
            print(feat.content)

if __name__ == "__main__":
    emain()

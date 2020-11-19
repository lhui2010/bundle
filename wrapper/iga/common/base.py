#!/usr/bin/env python

import argparse
import textwrap
import subprocess
from parse import *
import re
from collections import defaultdict
import base
import cfg 
import coloredlogs, logging

"""
Unfinished
"""

__all__ = ['JSONDecoder', 'JSONDecodeError']


# Create a logger object.
logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)
## Some examples.
#logger.debug("this is a debugging message")
#logger.info("this is an informational message")
#logger.warning("this is a warning message")
#logger.error("this is an error message")
#logger.critical("this is a critical message")

#def cmd_log_and_execute(cmd):
#    def prepare_workdir():
#        pass
#    logger.info(cmd)
#    subprocess.run(cmd, shell = True)


class DictDb():
    """
    The data structure for storing ctf files
    """
    def __init__(self):
        self.dictdb = defaultdict(dict)
        """
        whether this dictionary is nested (has sections) like 
        [abc]
        a=1
        b=2
        """
        self.has_section = False

    def append_val(self, key, val, section = ''):
        """
        Add tag, value to dict
        """
        if(section != ''):
            self.dictdb[section][key] = val
            self.has_section = True
        else:
            self.dictdb[key] = val

    def change_val(self, key, val, section = ''):
        """
        Add tag, value to dict
        """
        if(section != ''):
            self.dictdb[section][key] = val
        else:
            self.dictdb[key] = val

    def get_dict_text(self, seperator='='):
        """
        print this dict to text
        """
        return_text = ''
        if self.has_section:
            for section in self.dictdb:
                return_text += ('[{}]'.format(section) + "\n")
                for k in self.dictdb[section]:
                    return_text += ('{}{}{}'.format(k, seperator, self.dictdb[section][k]) + "\n")
        else:
            for k in self.dictdb:
                return_text += ('{}{}{}'.format(k, seperator, self.dictdb[k]) + "\n")
        return return_text


class Config():
    """
    Give arguments, will return the specific CTLs
    Support Direct Print
    Support tag value change
    """
    def __init__(self, cfg_type=""):
        self.content = ''
        self.dictdb = DictDb()
        self.seperator = ''
        if(cfg_type != ""):
            self.load(cfg_type)

    def load(self, cfg_type='falcon', seperator = '='):
        try:
            self.content = cfg.cfg[cfg_type]
            self.seperator = cfg.seperator[cfg_type]
        except KeyError as e:
            logger.error("Unknown type of cfg file: {}".format(cfg_type))

#Seperator for tag and value, like tag=value is default
        this_list = self.content.splitlines()
        for line in this_list:
            if(line.rstrip() == ''):
#Skip blank lines
                continue
            elif(line.startswith('[')):
#Finding section definition
                section = parse('[{}]', line)[0]
            elif(line.startswith('#')):
#Finding section annotation
                section_annotation = re.sub(r'^#+', '', line)
            else:
#Findng value assignments
                trimmed_line = re.sub('#.*', '', line)
                try:
                    (key, value) = trimmed_line.split(self.seperator)
                except ValueError as e:
                    logger.error("Split error on line: {}".format(line))
                (key, value) = (key.strip(), value.strip())
                if(section != ''):
                    self.dictdb.append_val(key=key, val=value, section=section)
                else:
                    self.dictdb.append_val(key, value)
            
    def update(self, args):
        """
        Different from change_val in dictdb, this allows input like :
        "[general]genomesize=12M;[general]threads=11"
        """
        mylist = args.split(';')
        for this_arg in mylist:
            section = ''
            if('[' in this_arg):
                (section, key, value) = parse('[{}]{}'+self.seperator+'{}', this_arg)
                self.dictdb.change_val(key=key, val=value, section=section)
            else:
                (key, value) = parse('{}'+self.seperator+'{}', this_arg)
                self.dictdb.change_val(key,value)

    def gettext(self):
        return self.dictdb.get_dict_text(self.seperator)

    @staticmethod
    def get_fofn(file_list, fofn_file):
        pass


def main():
    prog_name = "RunFalcon"
    usage = "Run Falcon With Fasta Input"

    parser = argparse.ArgumentParser(
        prog=prog_name, 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(usage), 
        epilog="")
    parser.add_argument("fasta", help="Raw assembly")
    args = parser.parse_args()  

    cfg = Config('falcon')
    cfg.update('[General]genome_size=10')
    cfg.update('[General]input_fofn=/dev/zero')
    logger.info(cfg.gettext())
    #falcon_run(args.fasta)
#    flanking_distance = args.flanking

if __name__ == "__main__":
    main()

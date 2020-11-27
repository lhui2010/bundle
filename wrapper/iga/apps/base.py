#!/usr/bin/env python

import argparse
import os
import shutil
import textwrap
import subprocess
from parse import *
import re
from collections import defaultdict
import iga.apps.cfg
import coloredlogs, logging
import os.path as op
import six
import sys

from iga.apps import cfg
from iga.utils.natsort import natsorted

"""
The basic library for iga, some of the functions are adapted from jcvi
"""

# Create a logger object.
logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)


## Some examples.
# logger.debug("this is a debugging message")
# logger.info("this is an informational message")
# logger.warning("this is a warning message")
# logger.error("this is an error message")
# logger.critical("this is a critical message")

def mkdir(dirname, overwrite=False):
    """
    Wraps around os.mkdir(), but checks for existence first.
    """
    if op.isdir(dirname):
        if overwrite:
            shutil.rmtree(dirname)
            os.mkdir(dirname)
            logging.debug("Overwrite folder `{0}`.".format(dirname))
        else:
            return False  # Nothing is changed
    else:
        try:
            os.mkdir(dirname)
        except:
            os.makedirs(dirname)
        logging.debug("`{0}` not found. Creating new.".format(dirname))

    return True


def split_fasta(fasta, workdir, chunk=100):
    file_list = sh('split_fastav3.pl {0} {2} && mv {0}._ {1}'.format(
        fasta, workdir, str(chunk)))
    return file_list

def sh(cmd, debug=False):
    """
    run command directly with subprocess.run
    :param cmd:
    :return:
    """
    ret = ''
    logger.info(cmd)
    prior_cmd = 'set -eo pipefail\n'
    if(debug == False):
        ret = subprocess.check_output(prior_cmd + cmd, shell=True).decode()
    return ret

def bsub(cmd, queue='Q104C512G_X4'):
    """
    submit jobs via bsub
    :param cmd:
    :return:
    """
    logger.info(cmd)
    bsub_cmd = 'bsub -q {}  -o output.%J -e error.%J '.format(queue)
    prior_cmd = 'set -eo pipefail\n'
    subprocess.run(bsub_cmd + '"' + prior_cmd + cmd + '"', shell=True)

conda_act = r"""
source ~/lh/anaconda3/etc/profile.d/conda.sh
conda activate {}
"""

workdir_sh = r"""
mkdir -p {0}
cd {0}
"""


def splitall(path):
    allparts = []
    while True:
        path, p1 = op.split(path)
        if not p1:
            break
        allparts.append(p1)
    allparts = allparts[::-1]
    return allparts


class ActionDispatcher(object):
    """
    jcvi
    This class will be invoked
    a) when the base package is run via __main__, listing all MODULESs
    a) when a directory is run via __main__, listing all SCRIPTs
    b) when a script is run directly, listing all ACTIONs
    This is controlled through the meta variable, which is automatically
    determined in get_meta().
    """

    def __init__(self, actions):

        self.actions = actions
        if not actions:
            actions = [(None, None)]
        self.valid_actions, self.action_helps = zip(*actions)

    def get_meta(self):
        args = sys.argv[0].split('/')[-3:]
        args[-1] = args[-1].replace(".py", "")
        if args[-2] == "iga":
            meta = "MODULE"
        elif args[-1] == "__main__":
            meta = "SCRIPT"
        else:
            meta = "ACTION"
        return meta, args

    def print_help(self):
        meta, args = self.get_meta()
        if meta == "MODULE":
            del args[0]
            args[-1] = meta
        elif meta == "SCRIPT":
            args[-1] = meta
        else:
            args[-1] += " " + meta

        help = "Usage:\n    python -m {0}\n\n\n".format(".".join(args))
        help += "Available {0}s:\n".format(meta)
        max_action_len = max(len(action) for action, ah in self.actions)
        for action, action_help in sorted(self.actions):
            action = action.rjust(max_action_len + 4)
            help += (
                    " | ".join((action, action_help[0].upper() + action_help[1:])) + "\n"
            )
        help += "\n"

        sys.stderr.write(help)
        sys.exit(1)

    def dispatch(self, globals):
        from difflib import get_close_matches

        meta = "ACTION"  # function is only invoked for listing ACTIONs
        if len(sys.argv) == 1:
            self.print_help()

        action = sys.argv[1]

        if not action in self.valid_actions:
            print("[error] {0} not a valid {1}\n".format(action, meta), file=sys.stderr)
            alt = get_close_matches(action, self.valid_actions)
            print(
                "Did you mean one of these?\n\t{0}\n".format(", ".join(alt)),
                file=sys.stderr,
            )
            self.print_help()

        globals[action](sys.argv[2:])


def dmain(mainfile, type="action"):
    cwd = op.dirname(mainfile)
    pyscripts = (
        [x for x in glob(op.join(cwd, "*", "__main__.py"))]
        if type == "module"
        else glob(op.join(cwd, "*.py"))
    )
    actions = []
    for ps in sorted(pyscripts):
        action = (
            op.basename(op.dirname(ps))
            if type == "module"
            else op.basename(ps).replace(".py", "")
        )
        if action[0] == "_":  # hidden namespace
            continue
        pd = get_module_docstring(ps)
        action_help = (
            [
                x.rstrip(":.,\n")
                for x in pd.splitlines(True)
                if len(x.strip()) > 10 and x[0] != "%"
            ][0]
            if pd
            else "no docstring found"
        )
        actions.append((action, action_help))

    a = ActionDispatcher(actions)
    a.print_help()


def get_module_docstring(filepath):
    "Get module-level docstring of Python module at filepath, e.g. 'path/to/file.py'."
    co = compile(open(filepath).read(), filepath, "exec")
    if co.co_consts and isinstance(co.co_consts[0], six.string_types):
        docstring = co.co_consts[0]
    else:
        docstring = None
    return docstring


def glob(pathname, pattern=None):
    """
    Wraps around glob.glob(), but return a sorted list.
    """
    import glob as gl

    if pattern:
        pathname = op.join(pathname, pattern)
    return natsorted(gl.glob(pathname))


class DictDb():
    """
    The data structure for storing ctf files
    Could be understand as a sub class of Config class
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

    def append_val(self, key, val, section=''):
        """
        Add tag, value to dict
        """
        if (section != ''):
            self.dictdb[section][key] = val
            self.has_section = True
        else:
            self.dictdb[key] = val

    def change_val(self, key, val, section=''):
        """
        Add tag, value to dict
        """
        if (section != ''):
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
        if (cfg_type != ""):
            self.load(cfg_type)

    def load(self, cfg_type='falcon', seperator='='):
        try:
            self.content = cfg.cfg[cfg_type]
            if cfg_type in cfg.seperator:
                self.seperator = cfg.seperator[cfg_type]
            else:
                self.seperator = seperator
        except KeyError as e:
            logger.error("Unknown type of cfg file: {}".format(cfg_type))

        # Seperator for tag and value, like tag=value is default
        this_list = self.content.splitlines()
        for line in this_list:
            if (line.rstrip() == ''):
                # Skip blank lines
                continue
            elif (line.startswith('[')):
                # Finding section definition
                section = parse('[{}]', line)[0]
            elif (line.startswith('#')):
                # Finding section annotation
                section_annotation = re.sub(r'^#+', '', line)
            else:
                # Findng value assignments
                trimmed_line = re.sub('#.*', '', line)
                try:
                    (key, value) = trimmed_line.split(self.seperator)
                except ValueError as e:
                    logger.error("Split error on line: {}".format(line))
                (key, value) = (key.strip(), value.strip())
                if (section != ''):
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
            if ('[' in this_arg):
                (section, key, value) = parse('[{}]{}' + self.seperator + '{}', this_arg)
                self.dictdb.change_val(key=key, val=value, section=section)
            else:
                (key, value) = parse('{}' + self.seperator + '{}', this_arg)
                self.dictdb.change_val(key, value)

    def get_text(self):
        return self.dictdb.get_dict_text(self.seperator)

    def write_to_file(self, output_file):
        with open(output_file, 'w') as fh:
            fh.write(self.dictdb.get_dict_text(self.seperator))

    @staticmethod
    def get_fofn(file_list, fofn_file):
        pass


def abspath_list(file_list):
    for i, v in enumerate(file_list):
        file_list[i] = op.abspath(v)

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
    # falcon_run(args.fasta)


#    flanking_distance = args.flanking

if __name__ == "__main__":
    main()

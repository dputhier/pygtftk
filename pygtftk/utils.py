"""A set of useful functions."""

import datetime
import glob
import io
import os
import random
import re
import string
import sys
import time
from collections import defaultdict, OrderedDict
from distutils.spawn import find_executable
from subprocess import PIPE
from subprocess import Popen
from tempfile import NamedTemporaryFile, mkdtemp

from pyparsing import Literal, CaselessLiteral, oneOf, nums, Word, Combine, Optional, operatorPrecedence, opAssoc, \
    Forward, ParseException

import pygtftk

# -------------------------------------------------------------------------
# VARIABLES
# -------------------------------------------------------------------------


# The level of verbosity
VERBOSITY = 0
WARN_REGION_SIZE = True
WARN_UNDEF = True

# A module wide list to keep track of temporary file.
TMP_FILE_LIST = []
TMP_DIR = None

# Temporary files
ADD_DATE = True

# ADD 'chr' when printing GTF
ADD_CHR = 0

# Whether chromosome file has been checked
CHROM_CHECKED = False

# Characters
TAB = '\t'
NEWLINE = '\n'

# Whether message should be written to a file
# Contains None or Argparse.Filetype('w')

MESSAGE_FILE = None


# ---------------------------------------------------------------
# Error class
# ---------------------------------------------------------------


class GTFtkError(Exception):

    def __init__(self, value):
        self.value = value

        if pygtftk.__NON_INTERACTIVE__:

            message(value, type="ERROR")
        else:
            raise GTFtkInteractiveError(value)

    def __str__(self):
        return repr(self.value)


class GTFtkInteractiveError(Exception):
    pass


# ---------------------------------------------------------------
# Directories
# ---------------------------------------------------------------

def mkdir_p(path):
    """Create a directory silently if already exists.

    :Example:

    >>> from pygtftk.utils import check_file_or_dir_exists
    >>> from pygtftk.utils import mkdir_p
    >>> mkdir_p('/tmp/test_gtftk_mkdir_p')
    >>> assert check_file_or_dir_exists('/tmp/test_gtftk_mkdir_p')
    """

    if not os.path.exists(path):
        if os.path.isfile(path):
            raise GTFtkError("The path provided is a file not a directory...")
        else:
            try:
                os.makedirs(path)
            except:
                raise GTFtkError("Unable to create directory : " + path)


def make_tmp_file(prefix='tmp',
                  suffix='',
                  store=True,
                  dir=None):
    """
    This function should be call to create a temporary file as all files
    declared in TMP_FILE_LIST will be remove upon command exit.

    :param prefix: a prefix for the temporary file (all of them will contain 'pygtftk').
    :param suffix: a suffix (default to '').
    :param store: declare the temporary file in utils.TMP_FILE_LIST. In pygtftk, \
    the deletion of these files upon exit is controled through -k.
    :param dir: a target directory.

    :Example:

    >>> from pygtftk.utils import make_tmp_file
    >>> tmp_file = make_tmp_file()
    >>> assert os.path.exists(tmp_file.name)
    >>> tmp_file = make_tmp_file(prefix="pref")
    >>> assert os.path.exists(tmp_file.name)
    >>> tmp_file = make_tmp_file(suffix="suf")
    >>> assert os.path.exists(tmp_file.name)

    """

    dir_target = None

    if dir is None:
        if TMP_DIR is not None:
            if not os.path.exists(TMP_DIR):
                msg = "Creating directory {d}."
                message(msg.format(d=TMP_DIR), type="INFO")
                os.mkdir(TMP_DIR)
                dir_target = TMP_DIR

    else:
        if not os.path.exists(dir):
            msg = "Creating directory {d}."
            message(msg.format(d=dir), type="INFO")
            os.mkdir(dir)
            dir_target = dir

    tmp_file = NamedTemporaryFile(delete=False,
                                  mode='w',
                                  prefix=prefix + "_pygtftk_",
                                  suffix=suffix,
                                  dir=dir_target)

    if store:
        TMP_FILE_LIST.append(tmp_file.name)

    return tmp_file


def make_tmp_dir(prefix='tmp',
                 store=True):
    """
    This function should be call to create a temporary file as all files
    declared in TMP_FILE_LIST will be remove upon command exit.

    :param prefix: a prefix for the temporary file (all of them will contain 'pygtftk').
    :param store: declare the temporary file in utils.TMP_FILE_LIST. In pygtftk, \
    the deletion of these files upon exit is controled through -k.

    :Example:

    >>> from pygtftk.utils import  make_tmp_dir
    >>> import os
    >>> import shutil
    >>> tp_dir = make_tmp_dir()
    >>> assert os.path.exists(tp_dir)
    >>> assert os.path.isdir(tp_dir)
    >>> shutil.rmtree(tp_dir)

    """

    dir_name = mkdtemp(prefix=prefix + "_pygtftk_")
    message("Creating directory : " + dir_name, type="DEBUG")

    if store:
        TMP_FILE_LIST.append(dir_name)

    return dir_name


# ---------------------------------------------------------------
# get the path to an example file
# ---------------------------------------------------------------

def get_example_feature():
    """
    Returns an example Feature (i.e GTF line).

    :Example:

    >>> from pygtftk.utils import get_example_feature
    >>> a = get_example_feature()
    """

    import pygtftk.Line

    alist = ['chr1', 'Unknown', 'transcript', 100, 200]
    alist += ['.', '+', '.', {'transcript_id': 'g1t1', 'gene_id': 'g1'}]

    return pygtftk.Line.Feature.from_list(alist)


def get_example_file(datasetname="simple", ext="gtf"):
    """
    Get the path to a gtf file example located in pygtftk library path.

    :param datasetname: The name of the dataset. One of 'simple', 'mini_real', 'simple_02, 'simple_03'.
    :param ext: Extension. For 'simple' dataset, can be one of 'bw', 'fa', 'fa.fai', 'chromInfo', 'bt*', 'fq', 'gtf' or '.*'.
    :Example:

    >>> from pygtftk.utils import get_example_file
    >>> a= get_example_file()
    >>> assert a[0].endswith('gtf')
    >>> a= get_example_file(ext="bam")
    >>> assert a[0].endswith('bam')
    >>> a= get_example_file(ext="bw")
    >>> assert a[0].endswith('bw')

    """

    """
    if datasetname == "mini_real" and ext == "bam":
        message("Retrieving BAM file")
        req = requests.get('https://tinyurl.com/y7ael4fo')
        with open("ENCFF742FDS_H3K4me3_K562_sub.bam", 'wb') as f:
            f.write(req.content)

        message("Retrieving BAI file")
        req = requests.get('https://tinyurl.com/y8j2rg2l')
        with open("ENCFF742FDS_H3K4me3_K562_sub.bam.bai", 'wb') as f:
            f.write(req.content)

        samfile = pysam.Samfile("ENCFF742FDS_H3K4me3_K562_sub.bam", "rb")

        return samfile
    """

    file_path = glob.glob(os.path.join(os.path.dirname(pygtftk.__file__),
                                       'data',
                                       datasetname,
                                       "*" + ext))

    return sorted(file_path)


# ---------------------------------------------------------------
# File processing
# ---------------------------------------------------------------


def add_prefix_to_file(infile, prefix=None):
    """ Returns a file object or path (as string) in which the
    base name is prefixed with 'prefix'. Returns a file object (write mode)
    or path depending on input type.

    :param infile: A file object or path.
    :param prefix: the prefix to add.

    :Example:

    >>> from pygtftk.utils import add_prefix_to_file
    >>> result = 'data/simple/bla_simple.gtf'
    >>> assert add_prefix_to_file('data/simple/simple.gtf', "bla_") == result

    """

    if prefix is None:
        return infile

    if isinstance(infile, io.IOBase):
        new_file = infile.name
    else:
        new_file = infile

    new_file = os.path.join(os.path.dirname(new_file),
                            prefix + os.path.basename(new_file))

    if isinstance(infile, io.IOBase):
        return open(new_file, "w")
    else:
        return new_file


def chomp(string):
    r"""
    Delete carriage return and line feed from end of string.

    :Example:

    >>> from pygtftk.utils import chomp
    >>> from pygtftk.utils import NEWLINE
    >>> assert NEWLINE not in chomp("blabla\r\n")

    """
    string = string.rstrip('\r\n')
    return string


def simple_line_count(afile):
    """Count the number of lines in a file.

    :param afile: A file object.

    :Example:

    >>> from pygtftk.utils import get_example_file
    >>> from pygtftk.utils import simple_line_count
    >>> my_file = get_example_file(datasetname="simple", ext="gtf")
    >>> my_file_h = open(my_file[0], "r")
    >>> assert simple_line_count(my_file_h) == 70

    """
    if afile.close:
        afile = open(afile.name, "r")

    lines = 0
    for _ in afile:
        lines += 1
    afile.close()

    return lines


def simple_nb_column(afile, separator="\t"):
    """Count the maximum number of columns in a file.

    :param separator: the separator (default "\t").
    :param afile: file name.

    :Example:

    >>> from pygtftk.utils import get_example_file
    >>> from pygtftk.utils import simple_nb_column
    >>> my_file = get_example_file(datasetname="simple", ext="gtf")
    >>> my_file_h = open(my_file[0], "r")
    >>> assert simple_nb_column(my_file_h) == 9

    """
    if afile.close:
        afile = open(afile.name, "r")

    nb_field = 0

    for line in afile:
        cur_nb_field = len(line.split(separator))
        if cur_nb_field > nb_field:
            nb_field = cur_nb_field

    afile.close()

    return nb_field


def head_file(afile=None, nb_line=6):
    """Display the first lines of a file (debug).

    :param afile: an input file.
    :param nb_line: the number of lines to print.

    :Example:

    >>> from pygtftk.utils import get_example_file
    >>> from pygtftk.utils import head_file
    >>> my_file = open(get_example_file()[0], "r")
    >>> #head_file(my_file)

    """

    if afile.close:
        afile = open(afile.name, "r")

    n = 1
    for line in afile:
        if n <= nb_line:
            print(chomp(line))
        n += 1

    afile.close()


def is_fasta_header(string):
    """Check if the line is a fasta header.

    :param string: a character string.

    :Example:

    >>> from pygtftk.utils import is_fasta_header
    >>> assert is_fasta_header(">DFTDFTD")

    """

    if string.startswith(">"):
        return True
    else:
        return False


def check_file_or_dir_exists(file_or_dir=None):
    """Check if a file/directory or a list of files/directories exist. Raise error if a file is not found.

    :param file_or_dir: file object or a list of file object.

    :Example:

    >>> from pygtftk.utils import get_example_file
    >>> from pygtftk.utils import check_file_or_dir_exists
    >>> assert check_file_or_dir_exists(get_example_file()[0])
    >>> assert check_file_or_dir_exists(get_example_file())

    """

    # Convert to a list
    if not isinstance(file_or_dir, list):
        file_or_dir = [file_or_dir]

    # Convert to filename
    file_or_dir = [x.name if isinstance(x, io.IOBase) else x for x in file_or_dir]

    for file_or_dir_cur in file_or_dir:

        if not os.path.exists(file_or_dir_cur):
            message("File not found: " + file_or_dir_cur, type="ERROR")

        else:
            message("Found file " + file_or_dir_cur, type="DEBUG")
    return True


def tab_line(token=None, newline=False):
    r"""Returns a tabulated string from an input list with an optional
    newline

    :param token: an input list.
    :param newline: if True a newline character is added to the string.

    :Example:

    >>> from pygtftk.utils import tab_line
    >>> from pygtftk.utils import TAB
    >>> assert tab_line(['a','b', 'c']) == 'a' + TAB + 'b' + TAB + 'c'

    """
    if not isinstance(token, list):
        raise GTFtkError("tab_line  needs a list as input.")

    token = [str(t) for t in token]

    if newline:
        return "\t".join(token) + "\n"

    return "\t".join(token)


def chrom_info_as_dict(chrom_info_file):
    """Check the format of the chromosome info file and return a dict.

    :param chrom_info_file: A file object.

    :Example:

    >>> from pygtftk.utils import get_example_file
    >>> from pygtftk.utils import chrom_info_as_dict
    >>> a = get_example_file(ext='chromInfo')[0]
    >>> b = chrom_info_as_dict(open(a, "r"))
    >>> assert b['chr1'] == 300
    >>> assert b['all_chrom'] == 900

    """

    message("Checking chromosome info file.",
            type="INFO")

    if chrom_info_file is None:
        raise GTFtkError("You must provide chromosome information file.")

    if chrom_info_file.closed:
        chrom_info_file = open(chrom_info_file.name, "r")

    chrom_len = OrderedDict()

    for line in chrom_info_file:

        if is_empty(line) or is_comment(
                line) or line.lower().startswith("chrom"):
            continue
        line = chomp(line)
        line = line.split("\t")

        try:
            chrom_len[line[0]] = int(line[1])

        except ValueError:
            continue

        if int(line[1]) <= 0:
            raise GTFtkError("Chromosome sizes should be greater than 0.")

    if len(list(chrom_len.keys())) == 0:
        raise GTFtkError("No chromosome length retrieved in chrom_info file. Is this file tabulated ?")

    genome_size = 0

    for _, value in list(chrom_len.items()):
        genome_size += value

    chrom_len["all_chrom"] = genome_size

    mes = "Found chromosome {i} (size: {size}) in {file})"

    if not pygtftk.utils.CHROM_CHECKED:

        for i in chrom_len:

            if i != "all_chrom":
                msg = mes.format(i=i,
                                 size=chrom_len[i],
                                 file=chrom_info_file.name)
                message(msg, type="DEBUG")
        pygtftk.utils.CHROM_CHECKED = True

    close_properly(chrom_info_file)

    return chrom_len


def chrom_info_to_bed_file(chrom_file, chr_list=None):
    """Return a bed file (file object) from a chrom info file.

    :param chrom_file: A file object.
    :param chr_list: A list of chromosome to be printed to bed file.

    >>> from  pygtftk.utils import chrom_info_to_bed_file
    >>> from  pygtftk.utils import get_example_file
    >>> from pybedtools import  BedTool
    >>> a = get_example_file(ext='chromInfo')
    >>> b = chrom_info_to_bed_file(open(a[0], 'r'))
    >>> c = BedTool(b.name)
    >>> chrs = []
    >>> for i in c: chrs += [i.chrom]
    >>> assert 'chr1' in chrs and 'chr2' in chrs
    >>> starts = []
    >>> for i in c: starts += [i.start]
    >>> assert 0 in starts
    >>> ends = []
    >>> for i in c: ends += [i.end]
    >>> assert 300 in ends and 600 in ends
    """

    message("Converting chrom info to bed format")

    out_file = make_tmp_file("chromInfo_", ".bed")
    chrom_dict = chrom_info_as_dict(chrom_file)

    for chrom, chrom_len in list(chrom_dict.items()):
        if chrom != "all_chrom":
            if chr_list is not None:
                if chrom in chr_list:
                    do_it = True
                else:
                    do_it = False
            else:
                do_it = True
        else:
            do_it = False

        if do_it:
            out_file.write("\t".join([chrom,
                                      "0",
                                      str(chrom_len),
                                      chrom,
                                      "0",
                                      "."]) + "\n")

    out_file.close()

    return open(out_file.name, "r")


def close_properly(*args):
    """Close a set of file if they are not None.

    :Example:

    >>> from pygtftk.utils import close_properly

    """
    for afile in args:
        if afile is not None:
            if afile != sys.stdout:
                afile.close()


def write_properly(a_string, afile):
    """Write a string to a file. If file is None, write string to stdout.

    :param a_string: a character string.
    :param afile: a file object.

    >>> from pygtftk.utils import write_properly

    """
    if afile is not None:

        if afile.name in ['<stdout>', '-']:
            sys.stdout.write(a_string + '\n')
        else:
            afile.write(a_string + '\n')


    else:
        sys.stdout.write(a_string + '\n')


def make_outdir_and_file(out_dir=None,
                         alist=None,
                         force=False):
    """Create output directory and a set of associated files (alist).
    Return a list of file handler (write mode) for the requested files.

    :param out_dir: the directory where file will be located.
    :param alist: the list of requested file names.
    :param force: if true, don't abort if directory already exists.

    :Example:

    >>> from pygtftk.utils import make_outdir_and_file

    """

    if out_dir is not None:
        # Working directory
        if os.path.exists(out_dir):
            if not force:
                raise GTFtkError("Aborting. Directory already exist.")
        else:
            os.makedirs(out_dir)

    out_fh_list = list()

    timestr = time.strftime("%Y%m%d-%H%M%S")

    for fn in alist:

        fn_split = os.path.splitext(fn)

        if pygtftk.utils.ADD_DATE:
            fn = fn_split[0] + "_" + timestr + fn_split[1]
        else:
            fn = fn_split[0] + fn_split[1]

        if out_dir is not None:
            out_fh_list.append(open(os.path.join(out_dir, fn), "w"))
        else:
            out_fh_list.append(open(fn, "w"))

    return out_fh_list


def silentremove(filename):
    """Remove silently (without error and warning) a file.

    :param filename: a file name (or file object).

    >>> from pygtftk.utils import silentremove
    >>> from pygtftk.utils import make_tmp_file
    >>> a = make_tmp_file()
    >>> assert os.path.exists(a.name)
    >>> silentremove(a.name)
    >>> assert not os.path.exists(a.name)

    """
    if isinstance(filename, io.IOBase):
        filename = filename.name

    try:
        os.remove(filename)
    except OSError as e:
        pass


def message(msg, nl=True, type="INFO", force=False):
    """Send a formated message on STDERR.

    :param msg: the message to be send.
    :param nl: if True, add a newline.
    :param type: Could be INFO, ERROR, WARNING.
    :param force: Force message, whatever the verbosity level.


    >>> from pygtftk.utils import message
    """

    now = datetime.datetime.now()

    if pygtftk.utils.VERBOSITY > 2:
        ho_min = str(now.hour) + ":" + str(now.minute).zfill(2) + ":" + str(now.second).zfill(2)
    else:
        ho_min = str(now.hour) + ":" + str(now.minute).zfill(2)

    do_it = False

    if type not in ["INFO", "ERROR", "WARNING", "DEBUG", "DEBUG_MEM"]:
        raise GTFtkError(
            "Type should be one of INFO, ERROR, WARNING, DEBUG, DEBUG_MEM.")

    if force:
        VERBOSITY_BACK = pygtftk.utils.VERBOSITY
        pygtftk.utils.VERBOSITY = 2

    if pygtftk.utils.VERBOSITY >= 2:
        do_it = True

    elif pygtftk.utils.VERBOSITY == 1:
        if type == "ERROR":
            do_it = True
        elif type == "WARNING":
            do_it = True
        elif type == "INFO":
            do_it = True
        elif type == "DEBUG_MEM":
            do_it = False
        elif type == "DEBUG":
            do_it = False

    elif pygtftk.utils.VERBOSITY == 0:
        if type == "ERROR":
            do_it = True
        elif type == "WARNING":
            do_it = True
        elif type == "INFO":
            do_it = False
        elif type == "DEBUG_MEM":
            do_it = False
        elif type == "DEBUG":
            do_it = False
    else:
        raise GTFtkError("Verbosity takes zero or positive integer values.")

    if do_it:
        msg = str(msg)

        if '__COMMAND__' in dir(pygtftk):
            cmd = pygtftk.__COMMAND__
        else:
            cmd = ""

        if cmd == "":
            if nl:
                msg = " |-- {c}-{a} : {b}\n".format(a=type,
                                                    b=msg,
                                                    c=ho_min)

            else:
                msg = " |-- {c}-{a} : {b}".format(a=type,
                                                  b=msg,
                                                  c=ho_min)

        else:
            if nl:
                msg = " |-- {c}-{a}-{d} : {b}\n".format(a=type,
                                                        b=msg,
                                                        c=ho_min,
                                                        d=cmd)
            else:
                msg = " |-- {c}-{a}-{d} : {b}".format(a=type,
                                                      b=msg,
                                                      c=ho_min,
                                                      d=cmd)
        if MESSAGE_FILE is None:
            sys.stderr.write(msg)
        else:
            MESSAGE_FILE.write(msg)
            MESSAGE_FILE.flush()

    if type == "ERROR":

        # Set sys.exit status to 1 to indicate an issue
        # to the system
        if pygtftk.__NON_INTERACTIVE__:
            message("Error encountered. System will exit after deleting temporary files.", type="DEBUG")

            if not pygtftk.__ARGS__['keep_all']:
                for i in flatten_list(TMP_FILE_LIST, outlist=[]):
                    message("ERROR encountered, deleting temporary file: " + i, type="DEBUG")
                    silentremove(i)
            else:
                message("Deletion of temporary files canceled by user.", type="DEBUG")
            sys.exit(1)
        else:
            raise GTFtkInteractiveError(msg)

    if force:
        pygtftk.utils.VERBOSITY = VERBOSITY_BACK


# ---------------------------------------------------------------
# Line type
# ---------------------------------------------------------------


def is_exon(string):
    """Does the string contains 'exon' preceded or not by spaces.

    :param string: The string to be tested.

    :Example:

    >>> from pygtftk.utils import is_exon
    >>> assert is_exon("Exon") and is_exon("exon")

    """
    if string.lower().strip() == "exon":
        return True

    return False


def is_empty(string):
    """
    Is the string empty.

    :param string: The string to be tested.

    :Example:

    >>> from pygtftk.utils import is_empty
    >>> assert is_empty("")

    """
    return not string.strip()


def is_comment(string):
    """Check wether a string is a comment (starts with #...).

    :param string: The string to be tested.

    :Example:

    >>> from pygtftk.utils import is_comment
    >>> assert is_comment("#bla")

    """

    return string.startswith("#")


# ---------------------------------------------------------------
# String
# ---------------------------------------------------------------

HEX_COLOR_REGEX = r'^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$'


def is_hex_color(input_string):
    regexp = re.compile(HEX_COLOR_REGEX)
    if regexp.search(input_string):
        return True
    return False


def random_string(n):
    """Returns a random string (alpha numeric) of length n."""
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(n))


def to_alphanum(string):
    """
    Returns a new string in which non alphanumeric character have been replaced by '_'. If non alphanumeric characters
    are found at the beginning or end of the string, they are deleted.

    :param string: A character string in which non alphanumeric char have to be replaced.


    :Example :

    >>> from pygtftk.utils import to_alphanum
    >>> assert to_alphanum("%gtf.bla") == 'gtf_bla'
    """

    if isinstance(string, int):
        string = str(string)

    if isinstance(string, list):
        string = "".join([str(x) for x in string])

    if string is None:
        string = 'NA'

    replacements = [('\W+', '_'), ('_+', '_')]

    for old, new in replacements:
        string = re.sub(old, new, string)

    string = string.rstrip("_")
    string = string.lstrip("_")
    return string


def left_strip_str(string):
    new_line = []

    for _, line in enumerate(string.split("\n")):
        if line == "":
            if not new_line:
                continue
        line = " " + re.sub("^\s+", "", line)
        new_line += [line.lstrip("\n")]

    return "\n".join(new_line)


# ---------------------------------------------------------------
# Lists and dicts
# ---------------------------------------------------------------

def to_list(obj, split_char=None):
    """ Convert a None, str, tuple to list. May also split if required."""
    if obj is None:
        obj = []
    elif isinstance(obj, str):
        if split_char is None:
            obj = [obj]
        else:
            obj = obj.split(split_char)
    elif isinstance(obj, tuple):
        obj = list(obj)
    return obj


def nested_dict(n, type):
    """"http://stackoverflow.com/questions/29348345"""
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n - 1, type))


def flatten_list(x, outlist=None):
    """Flatten a list of lists.

    :param x: a list or list of list.
    :param outlist: The output list.

    :Example:

    >>> from pygtftk.utils import flatten_list
    >>> b = ["A", "B", "C"]
    >>> a = ["a", "b", "c"]
    >>> c = ['a', 'b', 'c', 'A', 'B', 'C', 'string']
    >>> # Call it with an empty list (outlist)
    >>> # otherwise, element will be added to the existing
    >>> # outlist variable
    >>> assert flatten_list([a,b, "string"], outlist=[]) ==  c
    >>> assert flatten_list("string", outlist=[]) == ['string']

    """

    if outlist is None:
        outlist = []

    if not isinstance(x, (list, tuple)):
        outlist += [x]
    else:
        for i in x:
            outlist = flatten_list(i, outlist)
    return outlist


# See: https://stackoverflow.com/questions/716477/join-list-of-lists-in-python
def flatten_list_recur(a_list, sep=" "):
    """Join a list of list.
    :param a_list: a list.
    :param sep: the separator.

    """

    def _iter_flatten(root):

        if isinstance(root, (list, tuple)):
            for element in root:
                for e in _iter_flatten(element):
                    yield e
        else:
            yield root

    return sep.join(list(_iter_flatten(a_list)))


def sort_2_lists(list1, list2):
    """Sort list1 in increasing order and list2 occordingly.

    :param list1: the first list.
    :param list2: the second list.

    :Example:

    >>> from pygtftk.utils import sort_2_lists
    >>> a = ["c", "a", "b"]
    >>> b = ["A", "B", "C"]
    >>> c = sort_2_lists(a, b)
    >>> assert c == [('a', 'b', 'c'), ('B', 'C', 'A')]

    """

    tups = sorted(zip(list1, list2))

    return list(zip(*sorted(tups)))


# from http://stackoverflow.com/questions/8356501/python-format-tabular-output


def print_table(table):
    col_width = [max(len(x) for x in col) for col in zip(*table)]
    for line in table:
        print("| " + " | ".join("{:{}}".format(x, col_width[i])
                                for i, x in enumerate(line)) + " |")


def call_nested_dict_from_list(data, args=None):
    if args is None:
        args = []

    if args and data:
        element = args[0]
        if element:
            value = data.get(element)
        return value if len(args) == 1 else call_nested_dict_from_list(value, args[1:])


# ---------------------------------------------------------------
# Stats
# ---------------------------------------------------------------

def median_comp(alist):
    """Compute the median from a list.

    :param alist: a list.

    :Example:

    >>> from pygtftk.utils import median_comp
    >>> a = [10,20,30,40,50]
    >>> assert median_comp(a) == 30
    >>> a = [10,20,40,50]
    >>> assert median_comp(a) == 30

    """
    if len(alist) % 2 != 0:
        return sorted(alist)[len(alist) // 2]
    else:
        midavg = (sorted(alist)[len(alist) // 2] + sorted(alist)[len(alist) // 2 - 1]) / 2
        return midavg


def intervals(l, n, silent=False):
    """
    Returns a list of tuple that correspond to n ~ equally spaced intervals
    between 0 and l. Note that n must be lower than len(l).

    :param silent: use a silent mode.
    :param l: a range
    :param n: number of equally spaced intervals.

    :Example:

    >>> from pygtftk.utils import  intervals
    >>> l = range(10)
    >>> a = intervals(l, 3)
    >>> assert a == [(0, 3), (3, 6), (6, 10)]

    """

    if not l:
        raise GTFtkError("The list provided to 'intervals' function is empty.")

    if n <= 0 or n is None:
        raise GTFtkError("The n value provided to 'interval' function should be positive and not null.")

    result = list()

    if not len(l) >= n and not silent:
        raise GTFtkError("Cant' create " + str(n) + " equally spaced intervals "
                                                    "between 0 and " + str(l) + ".")

    def chunks(l, n):
        """ Yield n successive chunks from l.
        """

        newn = int(len(l) // n)
        for i in range(0, n - 1):
            yield l[i * newn:i * newn + newn + 1]
        yield l[n * newn - newn:]

    for i in chunks(l, n):
        result.append((i[0], i[len(i) - 1]))

    result[len(result) - 1] = (result[len(result) - 1][0],
                               result[len(result) - 1][1] + 1)
    return result


def check_r_installed():
    """Check R is installed.

    :Example:

    >>> from pygtftk.utils import check_r_installed
    >>> # check_r_installed()
    """
    if find_executable('R') is None:
        raise GTFtkError("R software was not found and is required.")


def check_r_packages(r_pkg_list=None, no_error=True):
    """
    Return True if R packages are installed. Return False otherwise.

    :param no_error: don't raise an error if some packages are not found.
    :param r_pkg_list: the list of R packages.

    :Example:

    >>> from pygtftk.utils import check_r_packages
    """

    p1 = Popen(["echo",
                "paste(rownames(installed.packages()), collapse=',')"],
               stdout=PIPE)
    p2 = Popen(["R", "--slave"], stdin=p1.stdout, stdout=PIPE)
    p1.stdout.close()
    installed, err = p2.communicate()
    installed = re.sub('^.*?"', '', installed)
    installed = re.sub('".*\n', '', installed)
    installed = installed.split(",")

    r_pkg_not_found = list(set(r_pkg_list) - set(installed))

    if len(r_pkg_not_found) > 0:
        message("Required R packages: " + " ".join(r_pkg_not_found),
                type="INFO")
        if no_error:
            message("The following R package(s) were not found on this system: " + ", ".join(r_pkg_not_found),
                    type="INFO")
        else:
            message("The following R package(s) were not found on this system: " + ", ".join(r_pkg_not_found),
                    type="ERROR")
        return False
    return True


# ---------------------------------------------------------------
# Check boolean expression
# ---------------------------------------------------------------

def check_boolean_exprs(exprs=None, operand=(), send_error=True):
    '''
    Check whether a boolean expression is properly formed.

    :param exprs: The string to evaluate.
    :param operand: The name of the operands.
    :param send_error: Whether to throw an error if expression is malformed.
    :return: A boolean.

    :Example:

    >>> from pygtftk.utils import check_boolean_exprs
    >>> assert check_boolean_exprs('s > 1 and (s < 2 or y < 2.5)', operand=['s', 'y'])

    '''
    lparen = Literal("(")
    rparen = Literal(")")
    and_operator = CaselessLiteral("and")
    or_operator = CaselessLiteral("or")
    comparison_operator = oneOf(['==', '!=', '>', '>=', '<', '<='])
    point = Literal('.')
    exponent = CaselessLiteral('E')
    plusorminus = Literal('+') | Literal('-')
    number = Word(nums)
    integer = Combine(Optional(plusorminus) + number)
    float_nb = Combine(integer +
                       Optional(point + Optional(number)) +
                       Optional(exponent + integer))
    value = float_nb
    identifier = oneOf(operand, caseless=False)  # .setParseAction(_embed)
    group_1 = identifier + comparison_operator + value
    group_2 = value + comparison_operator + identifier
    comparison = group_1 | group_2
    boolean_expr = operatorPrecedence(comparison,
                                      [(and_operator, 2, opAssoc.LEFT),
                                       (or_operator, 2, opAssoc.LEFT)])

    boolean_expr_par = lparen + boolean_expr + rparen

    expression = Forward()
    expression << boolean_expr | boolean_expr_par

    try:
        expression.parseString(exprs, parseAll=True)
        return True
    except ParseException as err:
        if send_error:
            message(err.msg, force=True)
            message('Operand should be one of: ' + ", ".join(operand))
            message("Boolean expression not supported.", type="ERROR")
        return False


# ---------------------------------------------------------------
# COLORS
# ---------------------------------------------------------------


ALL_MPL_PALETTES = ['viridis', 'plasma', 'inferno', 'magma',
                    'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
                    'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
                    'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
                    'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
                    'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
                    'hot', 'afmhot', 'gist_heat', 'copper',
                    'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
                    'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic',
                    'Pastel1', 'Pastel2', 'Paired', 'Accent',
                    'Dark2', 'Set1', 'Set2', 'Set3',
                    'tab10', 'tab20', 'tab20b', 'tab20c',
                    'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
                    'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',
                    'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar']

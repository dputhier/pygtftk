"""A set of useful functions."""

import datetime
import glob
import os
import re
import sys
import time
from collections import defaultdict
from distutils.spawn import find_executable
from subprocess import PIPE
from subprocess import Popen
from tempfile import NamedTemporaryFile

import pysam
import requests

import pygtftk

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

# R libraries
R_LIB = defaultdict(list)


# ---------------------------------------------------------------
# Directories
# ---------------------------------------------------------------

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


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

    if dir is None:
        if TMP_DIR is not None:
            if not os.path.exists(TMP_DIR):
                msg = "Creating directory {d}."
                message(msg.format(d=TMP_DIR), type="INFO")
                os.mkdir(TMP_DIR)

        tmp_file = NamedTemporaryFile(delete=False,
                                      prefix=prefix + "_gtftk_",
                                      suffix=suffix,
                                      dir=TMP_DIR)
    else:
        if not os.path.exists(dir):
            msg = "Creating directory {d}."
            message(msg.format(d=dir), type="INFO")
            os.mkdir(dir)

        tmp_file = NamedTemporaryFile(delete=False,
                                      prefix=prefix + "_gtftk_",
                                      suffix=suffix,
                                      dir=dir)

    if store:
        TMP_FILE_LIST.append(tmp_file.name)

    return tmp_file


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
    :param ext: Extension. For 'simple' dataset, can be one of 'bam', \
    'bam.bai', 'bw', 'fa', 'fa.fai', 'chromInfo', 'bt*', 'fq', 'gtf' or '.*'.
    :Example:

    >>> from pygtftk.utils import get_example_file
    >>> a= get_example_file()
    >>> assert a[0].endswith('gtf')
    >>> a= get_example_file(ext="bam")
    >>> assert a[0].endswith('bam')
    >>> a= get_example_file(ext="bw")
    >>> assert a[0].endswith('bw')

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

    file_path = glob.glob(os.path.join(os.path.dirname(pygtftk.__file__),
                                       'data',
                                       datasetname,
                                       "*" + ext))

    return file_path


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

    >>> from pygtftk.utils import get_example_file
    >>> from pygtftk.utils import add_prefix_to_file
    >>> gtf = get_example_file()[0]
    >>> result = ppygtftkk
    >>> assert add_prefix_to_file(get_example_file()[0], "bla_")[-32::1] == result

    """

    if prefix is None:
        return infile

    if isinstance(infile, file):
        new_file = infile.name
    else:
        new_file = infile

    new_file = os.path.join(os.path.dirname(new_file),
                            prefix + os.path.basename(new_file))

    if isinstance(infile, file):
        return open(new_file, "w")
    else:
        return new_file


def chomp(string):
    r"""
    Delete carriage return and line feed from end of string.

    :Example:

    >>> from pygtftk.utils import chomp
  pygtftk assert "\r" not in chomp("blabla\r\n")
    >>> assert "\n" not in chomp("blabla\r\n")

    """
    string = string.rstrip('\r\n')
    return string


def simple_line_count(afile):
    """Count the number of lines in a file.

    :param afile: A file object.

    :Example:

    >>> from pygtftk.utils import get_exampygtftkile
    >>> from pygtftk.utils import simple_lpygtftkount
    >>> my_file = get_example_file(datasetname="simple", ext="gtf")
    >>> my_file_h = open(my_file[0], "r")
    >>> assert simple_line_count(my_file_h) == 70

    """
    if afile.close:
        afile = open(afile.name, "r")

    lines = 0
    for line in afile:
        lines += 1
    afile.close()

    return lines


def simple_nb_column(afile, separator="\t"):
    """Count the maximum number of columns in a file.

    :param separator: the separator (default "\t").
    :param afile: file name.

    :Example:

    >>> from pygtftk.utils import get_exampygtftkile
    >>> from pygtftk.utils import simple_npygtftkumn
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

    >>> from pygtftk.utils import get_exampygtftkile
    >>> from pygtftk.utils import head_filpygtftk >>> my_file = open(get_example_file()[0], "r")
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

    >>> from pygtftk.utils import is_fastapygtftker
    >>> assert is_fasta_header(">DFTDFTD")

    """

    if string.startswith(">"):
        return True
    else:
        return False


def check_file_or_dir_exists(file_or_dir=None):
    """Check if a file/directory or a list of files/directories exist. Raise error if a file is not found.

    :param file_or_dir: file object or a list of file object.
    :param  is_list: Is file_or_dir a list ?

    :Example:

    >>> from pygtftk.utils import get_exampygtftkile
    >>> from pygtftk.utils import check_fipygtftk_dir_exists
    >>> assert check_file_or_dir_exists(get_example_file()[0])
    >>> assert check_file_or_dir_exists(get_example_file())

    """

    test_list = []

    # Convert to a list
    if not isinstance(file_or_dir, list):
        file_or_dir = [file_or_dir]

    # Convert to filename
    file_or_dir = [x.name if isinstance(x, file) else x for x in file_or_dir]

    for file_or_dir_cur in file_or_dir:
        if not os.path.exists(file_or_dir_cur):
            raise pygtftk.error.GTFtkError("File not found: " + file_or_dir_cur)

        else:
            message("Found file " + file_or_dir_cur)
    return True


def tab_line(token=None, newline=False):
    r"""Returns a tabulated string from an input list with an optional
    newline

    :param token: an input list.
    :param newline: if True a newline character is added to the string.

    :Example:

    >>> from pygtftk.utils import tab_linepygtftk>>> assert tab_line(["a","b", "c"]) == "a\tb\tc"

    """
    if not isinstance(token, list):
        raise pygtftk.error.GTFtkError("tab_line  eeds a list as input.")

    token = [str(t) for t in token]

    if newline:
        return "\t".join(token) + "\n"

    return "\t".join(token)


def chrom_info_as_dict(chrom_info_file):
    """Check the format of the chromosome info file and return a dict.

    :param chrom_info_file: A file object.

    :Example:

    >>> from pygtftk.utils import get_exampygtftkile
    >>> from pygtftk.utils import chrom_inpygtftk_dict
    >>> a = get_example_file(ext='chromInfo')[0]
    >>> b = chrom_info_as_dict(open(a, "r"))
    >>> assert b['chr1'] == 300
    >>> assert b['all_chrom'] == 900

    """

    message("Checking chromosome info file.",
            type="INFO")

    if chrom_info_file is None:
        raise pygtftk.error.GTFtkError("You must provide chromosome information file.")

    if chrom_info_file.closed:
        chrom_info_file = open(chrom_info_file.name, "r")

    chrom_len = defaultdict(int)

    for line in chrom_info_file:

        if is_empty(line) or is_comment(
                line) or line.lower().startswith("chrom"):
            continue
        line = chomp(line)
        line = line.split("\t")

        try:
            chrom_len[line[0]] = int(line[1])

        except:
            continue

        if int(line[1]) <= 0:
            raise pygtftk.error.GTFtkError("Chromosome sizes should be greater than 0.")

    if len(chrom_len.keys()) == 0:
        raise pygtftk.error.GTFtkError("No chromosome length retrieved in chrom_info file. Is this file tabulated ?")

    genome_size = 0

    for chrom, value in chrom_len.items():
        genome_size += value

    chrom_len["all_chrom"] = genome_size

    mes = "Found chromosome {i} (size: {size}) in {file})"

    if not pygtftk.utils.CHROM_CHECKED:

        for i in chrom_len:

            if i != "all_chrom":
                msg = mes.format(i=i,
                                 size=chrom_len[i],
                                 file=chrom_info_file.name)
                message(msg,
                        type="INFO")
        pygtftk.utils.CHROM_CHECKED = True

    close_properly(chrom_info_file)

    return chrom_len


def chrom_info_to_bed_file(chrom_file, chr_list=None):
    """Return a bed file (file object) from a chrom info file.

    :param chrom_file: A file object.
    :param chr_list: A list of chromosome to be printed to bed file.

    >>> from  pygtftk.utils import chrom_inpygtftk_bed_file
    >>> from  pygtftk.utils import get_exampygtftkile
    >>> from pybedtools import  BedTool
    >>> a = get_example_file(ext='chromInfo')
    >>> b = chrom_info_to_bed_file(open(a[0], 'r'))
    >>> c = BedTool(b.name)
    >>> d = c.__iter__()
    >>> i = d.next()
    >>> assert i.start == 0
    >>> assert i.end == 600
    >>> assert (i.end - i.start) == 600
    >>> i = d.next()
    >>> assert (i.end - i.start) == 300

    """

    message("Converting chrom info to bed format")

    out_file = make_tmp_file("chromInfo_", ".bed")
    chrom_dict = chrom_info_as_dict(chrom_file)

    for chrom, chrom_len in chrom_dict.items():
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

    >>> from pygtftk.utils import close_prpygtftky

    """
    for afile in args:
        if afile is not None:
            if afile != sys.stdout:
                afile.close()


def write_properly(string, afile):
    """Write a string to a file. If file is None, write string to stdout.

    :param string: a character string.
    :param afile: a file object.

    >>> from pygtftk.utils import write_prpygtftky

    """
    if afile is not None and afile.name != '<stdout>':
        afile.write(string + "\n")
    else:
        sys.stdout.write(string + "\n")


def make_outdir_and_file(out_dir=None,
                         alist=None,
                         force=False):
    """Create output directory and a set of associated files (alist).
    Return a list of file handler (write mode) for the requested files.

    :param out_dir: the directory where file will be located.
    :param alist: the list of requested file names.
    :param force: if true, don't abort if directory already exists.

    :Example:

    >>> from pygtftk.utils import make_outpygtftknd_file

    """

    if out_dir is not None:
        # Working directory
        if os.path.exists(out_dir):
            if not force:
                raise pygtftk.error.GTFtkError("Aborting. Directory already exist.")
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

    >>> from pygtftk.utils import silentrepygtftk    >>> from pygtftk.utils import make_tmppygtftk
    >>> a = make_tmp_file()


    """
    if isinstance(filename, file):
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
    ho_min = str(now.hour) + ":" + str(now.minute).zfill(2)
    do_it = False

    if type not in ["INFO", "ERROR", "WARNING", "DEBUG", "DEBUG_MEM"]:
        raise pygtftk.error.GTFtkError(
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
        raise pygtftk.error.GTFtkError("Verbosity takes zero or positive integer values.")

    if do_it:
        msg = str(msg)

        if '__COMMAND__' in dir(pygtftk):
            cmd = pygtftk.__COMMAND__
        else:
            cmd = ""

        if cmd == "":
            if nl:
                sys.stderr.write("    |--- {c}-{a} : {b}\n".format(a=type,
                                                                   b=msg,
                                                                   c=ho_min))
            else:
                sys.stderr.write("    |--- {c}-{a} : {b}".format(a=type,
                                                                 b=msg,
                                                                 c=ho_min))
        else:
            if nl:
                sys.stderr.write("    |--- {c}-{a}-{d} : {b}\n".format(a=type,
                                                                       b=msg,
                                                                       c=ho_min,
                                                                       d=cmd))
            else:
                sys.stderr.write("    |--- {c}-{a}-{d} : {b}".format(a=type,
                                                                     b=msg,
                                                                     c=ho_min,
                                                                     d=cmd))

    if type == "ERROR":

        # Set sys.exit status to 1 to indicate an issue
        # to the system

        message("Error encountered. System will exit after deleting temporary files.", type="DEBUG")

        if not pygtftk.__ARGS__['keep_all']:
            for i in flatten_list(TMP_FILE_LIST):
                message("ERROR encountered, deleting temporary file: " + i, type="DEBUG")
                silentremove(i)
        else:
            message("Deletion of temporary files canceled by user.", type="DEBUG")
        sys.exit(1)

    if force:
        pygtftk.utils.VERBOSITY = VERBOSITY_BACK


# ---------------------------------------------------------------
# Line type
# ---------------------------------------------------------------

def add_r_lib(cmd=None, libs=None):
    """Declare a new set of required R libraries.

    :param libs: Comma separated list of R libraries.
    :param cmd: The target command.

    """

    for lib in libs.split(","):
        R_LIB[cmd] += [lib]


# ---------------------------------------------------------------
# Line type
# ---------------------------------------------------------------


def is_exon(string):
    """Does the string contains 'exon' preceded or not by spaces.

    :param string: The string to be tested.

    :Example:

    >>> from pygtftk.utils import is_exon
pygtftk>> assert is_exon("Exon") and is_exon("exon")

    """
    if string.lower().strip() == "exon":
        return True

    return False


def is_empty(string):
    """
    Is the string empty.

    :param string: The string to be tested.

    :Example:

    >>> from pygtftk.utils import is_emptypygtftk>>> assert is_empty("")

    """
    return not string.strip()


def is_comment(string):
    """Check wether a string is a comment (starts with #...).

    :param string: The string to be tested.

    :Example:

    >>> from pygtftk.utils import is_commepygtftk  >>> assert is_comment("#bla")

    """

    return string.startswith("#")


# ---------------------------------------------------------------
# String
# ---------------------------------------------------------------

def to_alphanum(string):
    """
    Returns a new string in which non alphanumeric character have been replaced by '_'.

    :param string: A character string in which non alphanumeric char have to be replaced.
    :param replacement_list:
    :return:
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

    for pos, line in enumerate(string.split("\n")):
        if line == "":
            if new_line == []:
                continue
        line = " " + re.sub("^\s+", "", line)
        new_line += [line.lstrip("\n")]

    return "\n".join(new_line)


# ---------------------------------------------------------------
# Lists and dicts
# ---------------------------------------------------------------


def flatten_list(alist):
    """Flatten a list of lists.

    :param alist: a list or list of list.

    :Example:

    >>> b = ["A", "B", "C"]
    >>> a = ["a", "b", "c"]
    >>> c = ['a', 'b', 'c', 'A', 'B', 'C', 'string']
    >>> assert flatten_list([a,b, "string"]) ==  c
    >>> assert flatten_list("string") == ['string']

    """

    import itertools

    if isinstance(alist, list):
        return list(itertools.chain.from_iterable(
            itertools.repeat(x, 1) if isinstance(x, str) else x for x in alist))
    elif isinstance(alist, str):
        return [alist]
    else:
        raise pygtftk.error.GTFtkError("Should be a list or str.")


# See: https://stackoverflow.com/questions/716477/join-list-of-lists-in-python
def flatten_list_recur(a_list, sep=" "):
    """Join a list of list.
    :param a_list: a list.
    :param sep: the separator.

    """

    def _iterFlatten(root):

        if isinstance(root, (list, tuple)):
            for element in root:
                for e in _iterFlatten(element):
                    yield e
        else:
            yield root

    return sep.join(list(_iterFlatten(a_list)))


def sort_2_lists(list1, list2):
    """Sort list1 in increasing order and list2 occordingly.

    :param list1: the first list.
    :param list2: the second list.

    :Example:

    >>> from pygtftk.utils import sort_2_lpygtftk    >>> a = ["c", "a", "b"]
    >>> b = ["A", "B", "C"]
    >>> c = sort_2_lists(a, b)
    >>> assert c == [('a', 'b', 'c'), ('B', 'C', 'A')]

    """

    tups = sorted(zip(list1, list2))
    return zip(*sorted(tups))


# from http://stackoverflow.com/questions/8356501/python-format-tabular-output


def print_table(table):
    col_width = [max(len(x) for x in col) for col in zip(*table)]
    for line in table:
        print "| " + " | ".join("{:{}}".format(x, col_width[i])
                                for i, x in enumerate(line)) + " |"


def call_nested_dict_from_list(data, args=[]):
    if args and data:
        element = args[0]
        if element:
            value = data.get(element)
        return value if len(args) == 1 else call_nested_dict_from_list(value, args[1:])


# ---------------------------------------------------------------
# Stats
# ---------------------------------------------------------------

def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        https://stackoverflow.com/questions/8930370/where-can-i-find-mad-mean-absolute-deviation-in-scipy
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))

def median_comp(alist):
    """Compute the median from a list.

    :param alist: a list.

    :Example:

    >>> from pygtftk.utils import median_cpygtftk   >>> a = [10,20,30,40,50]
    >>> assert median_comp(a) == 30
    >>> a = [10,20,40,50]
    >>> assert median_comp(a) == 30

    """
    if len(alist) % 2 != 0:
        return sorted(alist)[len(alist) / 2]
    else:
        midavg = (sorted(alist)[len(alist) / 2] + sorted(alist)[len(alist) / 2 - 1]) / 2.0
        return midavg


def intervals(l, n, silent=False):
    """
    Returns a list of tuple that correspond to n ~ equally spaced intervals
    between 0 and l. Note that n must be lower than len(l).

    :param silent: use a silent mode.
    :param l: a range
    :param n: number of equally spaced intervals.

    :Example:

    >>> from pygtftk.utils import  intervapygtftk  >>> l = range(10)
    >>> a = intervals(l, 3)
    >>> assert a == [(0, 3), (3, 6), (6, 10)]

    """

    if not l:
        raise pygtftk.error.GTFtkError("The list provided to 'intervals' function is empty.")

    if n <= 0 or n is None:
        raise pygtftk.error.GTFtkError("The n value provided to 'interval' function should be positive and not null.")

    result = list()

    if not len(l) >= n and not silent:
        raise pygtftk.error.GTFtkError("Cant' create " + str(n) + " equally spaced intervals "
                                                                  "between 0 and " + str(l) + ".")

    def chunks(l, n):
        """ Yield n successive chunks from l.
        """

        newn = int(len(l) / n)
        for i in xrange(0, n - 1):
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

    >>> from pygtftk.utils import check_r_pygtftklled
    """
    if find_executable('R') is None:
        raise pygtftk.error.GTFtkError("R software was not found and is required.")


def check_r_packages(r_pkg_list=None, no_error=True):
    """
    Return True if R packages are installed. Return False otherwise.

    :param no_error: don't raise an error if some packages are not found.
    :param r_pkg_list: the list of R packages.

    :Example:

    >>> from pygtftk.utils import check_r_pygtftkges
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

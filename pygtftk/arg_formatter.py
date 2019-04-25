# -*- coding: utf-8 -*-
"""
Command Line Interface display format.
"""

import argparse
import glob
import io
import operator
import os
import re
import sys

import pybedtools
from pybedtools import BedTool

import pygtftk
from pygtftk.utils import check_file_or_dir_exists, make_tmp_file, close_properly
from pygtftk.utils import chrom_info_as_dict
from pygtftk.utils import message


# ---------------------------------------------------------------
# ArgFormatter class
# ---------------------------------------------------------------


class ArgFormatter(argparse.HelpFormatter):
    """
    A correction to the argument formatter. This ensure proper width
    of the first column when print_usage() is called.
    """

    def __init__(self, prog):
        super(ArgFormatter, self).__init__(prog,
                                           indent_increment=1,
                                           max_help_position=40,
                                           width=300)
        self._width = 1000

    def add_argument(self, action):
        if action.help is not argparse.SUPPRESS:

            # find all invocations
            get_invocation = self._format_action_invocation
            invocations = [get_invocation(action)]
            current_indent = self._current_indent
            for subaction in self._iter_indented_subactions(action):
                # compensate for the indent that will be added
                indent_chg = self._current_indent - current_indent
                added_indent = 'x' * indent_chg
                invocations.append(added_indent + get_invocation(subaction))
            # print('inv', invocations)

            # update the maximum item length
            invocation_length = max([len(s) for s in invocations])
            action_length = invocation_length + self._current_indent
            self._action_max_length = max(self._action_max_length,
                                          action_length)

            # add the item to the list
            self._add_item(self._format_action, [action])

    def _get_help_string(self, action):
        help_str = action.help
        if '%(default)' not in action.help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    if isinstance(action.default, io.IOBase):
                        help_str += ' (default: ' + \
                                    str(action.default.name) + ')'
                    else:
                        help_str += ' (default: %(default)s)'
        return help_str

    def _fill_text(self, text, width, indent):
        return text

    def _get_default_metavar_for_optional(self, action):
        return action.dest.lower()

    def _get_default_metavar_for_positional(self, action):
        return action.dest.lower()

    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            metavar, = self._metavar_formatter(action, default)(1)
            return metavar

        else:
            parts = []

            # if the Optional doesn't take a value, format is:
            #    -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)

            # if the Optional takes a value, format is:
            #    -s ARGS, --long ARGS
            else:
                default = self._get_default_metavar_for_optional(action)
                for option_string in action.option_strings:
                    parts.append(option_string)

                return ', '.join(parts)
            return ', '.join(parts)

    def _format_usage(self, usage, actions, groups, prefix):
        """From argparse.py"""
        if prefix is None:
            prefix = '\n  Usage: '

        # if usage is specified, use that
        if usage is not None:
            usage = usage % dict(prog=self._prog)

        # if no optionals or positionals are available, usage is just prog
        elif usage is None and not actions:
            usage = '%(prog)s' % dict(prog=self._prog)

        # if optionals and positionals are available, calculate usage
        elif usage is None:
            prog = '%(prog)s' % dict(prog=self._prog)

            # split optionals from positionals
            optionals = []
            positionals = []
            for action in actions:
                if action.option_strings:
                    optionals.append(action)
                else:
                    positionals.append(action)

            # build full usage string
            format_fun = self._format_actions_usage
            action_usage = format_fun(optionals + positionals, groups)
            usage = ' '.join([s for s in [prog, action_usage] if s])

            # wrap the usage parts if it's too long
            text_width = self._width - self._current_indent
            if len(prefix) + len(usage) > text_width:

                # break usage into wrappable parts
                part_regexp = r'\(.*?\)+|\[.*?\]+|\S+'
                opt_usage = format_fun(optionals, groups)
                pos_usage = format_fun(positionals, groups)
                opt_parts = re.findall(part_regexp, opt_usage)
                pos_parts = re.findall(part_regexp, pos_usage)

                # helper for wrapping lines
                def get_lines(parts, indent, prefix=None):
                    lines = []
                    line = []
                    if prefix is not None:
                        line_len = len(prefix) - 1
                    else:
                        line_len = len(indent) - 1
                    for part in parts:
                        if line_len + 1 + len(part) > text_width:
                            lines.append(indent + ' '.join(line))
                            line = []
                            line_len = len(indent) - 1
                        line.append(part)
                        line_len += len(part) + 1
                    if line:
                        lines.append(indent + ' '.join(line))
                    if prefix is not None:
                        lines[0] = lines[0][len(indent):]
                    return lines

                # if prog is short, follow it with optionals or positionals
                if len(prefix) + len(prog) <= 0.75 * text_width:
                    indent = ' ' * (len(prefix) + len(prog) + 1)
                    if opt_parts:
                        lines = get_lines([prog] + opt_parts, indent, prefix)
                        lines.extend(get_lines(pos_parts, indent))
                    elif pos_parts:
                        lines = get_lines([prog] + pos_parts, indent, prefix)
                    else:
                        lines = [prog]

                # if prog is long, put it on its own line
                else:
                    indent = ' ' * len(prefix)
                    parts = opt_parts + pos_parts
                    lines = get_lines(parts, indent)
                    if len(lines) > 1:
                        lines = []
                        lines.extend(get_lines(opt_parts, indent))
                        lines.extend(get_lines(pos_parts, indent))
                    lines = [prog] + lines

                # join lines into usage
                usage = '\n'.join(lines)

        # prefix with 'usage:'
        return '%s%s\n\n' % (prefix, usage)


# ---------------------------------------------------------------
# Ensure a numeric is in range
# ---------------------------------------------------------------

def get_truth(inp, relate, cut):
    ops = {'>': operator.gt,
           '<': operator.lt,
           '>=': operator.ge,
           '<=': operator.le,
           '=': operator.eq}
    return ops[relate](inp, cut)


def ranged_num(lowest=-1, highest=1, val_type=("int", "float"),
               linc=False, hinc=False):
    """Check a numeric is in expected range.

    :param lowest: The lowest accepted value.
    :param highest: The highest accepted value.
    :param val_type: A float of int.
    :param linc: The lowest value is included (<=).
    :param hinc: The highest value is included (>=).

    """

    if linc:
        lt = '<'
    else:
        lt = '<='

    if hinc:
        gt = '>'
    else:
        gt = '>='

    if isinstance(val_type, tuple):
        val_type = val_type[0]

    def type_func(a_value):
        if val_type == "int":
            try:
                a_value = int(a_value)
            except:
                raise argparse.ArgumentTypeError("Not an int.")
        else:
            try:
                a_value = float(a_value)
            except ValueError:
                raise argparse.ArgumentTypeError("Not a float.")

        if lowest is None:
            if highest is None:
                return a_value
            else:
                if get_truth(a_value, gt, highest):
                    raise argparse.ArgumentTypeError("Value not in range.")
        else:
            if highest is None:
                if get_truth(a_value, lt, lowest):
                    raise argparse.ArgumentTypeError("Value not in range.")
            else:
                if get_truth(a_value, lt, lowest) or get_truth(a_value, gt, highest):
                    raise argparse.ArgumentTypeError("Value not in range.")

        return a_value

    return type_func


# ---------------------------------------------------------------
# Check the chrom file format
# ---------------------------------------------------------------

class CheckChromFile(argparse.Action):
    """
    Check the chromosome file exists and has the proper format.
    """

    def __init__(self,
                 option_strings,
                 dest,
                 nargs=None,
                 const=None,
                 default=None,
                 type=None,
                 choices=None,
                 required=False,
                 help=None,
                 metavar=None):
        argparse.Action.__init__(self,
                                 option_strings=option_strings,
                                 dest=dest,
                                 nargs=nargs,
                                 const=const,
                                 default=default,
                                 type=type,
                                 choices=choices,
                                 required=required,
                                 help=help,
                                 metavar=metavar,
                                 )

    def __call__(self,
                 parser,
                 namespace,
                 values,
                 option_string=None):
        if values in ["mm8", "mm9", "mm10", "hg19", "hg38", "rn3", "rn4"]:
            chr_size = pybedtools.helpers.chromsizes(values)
            ## Delete haplotype chromosome
            ## unplaced contig and unlocalized contig
            regexp = re.compile('(_random)|(^chrUn)|(_hap\d+)|(_alt)|(^chrM$)')
            chr_size = {key: chr_size[key] for key in chr_size if not regexp.search(key)}
            tmp_file_chr = make_tmp_file(prefix='chromsize', suffix='.txt')
            for chrom, size in chr_size.items():
                tmp_file_chr.write(chrom + "\t" + str(size[1]) + "\n")
            tmp_file_chr.close()
            values = open(tmp_file_chr.name, 'r')

        else:
            check_file_or_dir_exists(values)
            values = open(values, "r")
            chrom_info_as_dict(values)

        # Add the attribute
        setattr(namespace, self.dest, values)


# ---------------------------------------------------------------
# Check file extension and format
# ---------------------------------------------------------------

class FormattedFile(argparse.FileType):
    """
    Check file extensions and format.

    :param mode: the mode ('r'...).
    :param mode: A string or tuple, The accepted file_ext  ('bed', 'bed.gz', 'txt', 'txt.gz', 'gtf', 'gtf.gz', 'fasta',
    'fasta.gz', 'zip', 'bigwig')

    """

    def __init__(self, mode='r', file_ext='bed', **kwargs):
        super(FormattedFile, self).__init__(mode, **kwargs)
        self.file_ext = file_ext

    def __call__(self, string):

        # ---------------------------------------------------------------
        # Check file extension
        # ---------------------------------------------------------------

        fasta_format_1 = '(\.[Ff][Aa][Ss][Tt][Aa]$)|(\.[Ff][Nn][Aa]$)'
        fasta_format_2 = '|(\.[Ff][Aa]$)|(\.[Ff][Aa][Ss]$)|(\.[Ff][Ff][Nn]$)|(\.[Ff][Rr][Nn]$)'
        fasta_regexp = fasta_format_1 + fasta_format_2
        fasta_regexp_gz = re.sub("\$", "\.[Gg][Zz]$", fasta_regexp)
        bed_regexp = '\.[Bb][Ee][Dd][3456]{0,1}$'
        bed_regexp_gz = re.sub("\$", "\.[Gg][Zz]$", bed_regexp)
        gtf_regexp = '\.[Gg][Tt][Ff]$'
        gtf_regexp_gz = re.sub("\$", "\.[Gg][Zz]$", gtf_regexp)
        txt_regexp = '(\.[Tt][Xx][Tt]$)|(\.[Cc][Ss][Vv]$)|(\.[Dd][Ss][Vv]$)|(\.[Tt][Aa][Bb]$)|(\.[Tt][Ss][Vv]$)'
        txt_regexp_gz = re.sub("\$", "\.[Gg][Zz]$", txt_regexp)
        bigwig_regexp = '(\.[Bb][Ww]$)|(\.[Bb][Ii][Gg][Ww][Ii][Gg]$)'
        zip_regexp = '\.[Zz][Ii][Pp]$'
        pdf_regexp = '\.[Pp][Dd][Ff]$'

        ext2regexp = {'bed': bed_regexp,
                      'bed.gz': bed_regexp_gz,
                      'gtf': gtf_regexp,
                      'gtf.gz': gtf_regexp_gz,
                      'fasta': fasta_regexp,
                      'fasta.gz': fasta_regexp_gz,
                      'txt': txt_regexp,
                      'txt.gz': txt_regexp_gz,
                      'bigwig': bigwig_regexp,
                      'zip': zip_regexp,
                      'pdf': pdf_regexp}

        # Set verbosity system wide as depending on
        # command line argument order, VERBOSITY (-V) can
        # be evaluated later...
        if '-V' in sys.argv:
            sys_args = ' '.join(sys.argv)
            verbosity_val = re.search('-V ?([01234])?', sys_args)
            if verbosity_val:
                pygtftk.utils.VERBOSITY = int(verbosity_val.group(1))
            else:
                pygtftk.utils.VERBOSITY = 0

        match = False

        if isinstance(self.file_ext, str):
            extension_list = [self.file_ext]
        else:
            extension_list = list(self.file_ext)

        for this_ext in extension_list:
            if re.search(ext2regexp[this_ext], string):
                match = True
                break

        if not match:
            message('Not a valid filename extension :' + string, type="WARNING")
            message('Extension expected: ' + ext2regexp[this_ext], type="ERROR")
            sys.exit()

        # ---------------------------------------------------------------
        # Check directory
        # ---------------------------------------------------------------

        outputdir = os.path.dirname(os.path.abspath(string))

        if not os.path.exists(outputdir):
            if 'w' in self._mode:
                message("Directory not found. Creating.", type="WARNING")
                os.makedirs(outputdir)

        # ---------------------------------------------------------------
        # Check format
        # ---------------------------------------------------------------

        # if bed3, bed4, bad5 convert to bed6

        if self._mode == 'r':
            if self.file_ext == 'bed':

                message("Checking BED file format (" + string + ").",
                        type="INFO")

                try:
                    file_bo = BedTool(string)
                    nb_line = len(file_bo)
                except:
                    msg = "Unable to load file: " + string + "."
                    message(msg, type="ERROR")
                    sys.exit()

                if nb_line == 0:
                    msg = "It seems that file " + string + " is empty."
                    message(msg, type="ERROR")
                    sys.exit()

                if file_bo.file_type != 'bed':
                    msg = "File {f} is not a valid bed file."
                    msg = msg.format(f=string)
                    message(msg, type="ERROR")
                    sys.exit()

                region_nb = 0
                field_count = file_bo.field_count()

                if field_count != 6:
                    message("Converting to bed6 format (" + string + ").", type="WARNING")
                    tmp_file = make_tmp_file(prefix="bed6_",
                                             suffix=".bed")
                    for record in file_bo:
                        region_nb += 1

                        if field_count < 4:
                            name = 'region_' + str(region_nb)
                        else:
                            name = record.name

                        fields = record.fields[0:3]
                        fields += [name, '0', '.']

                        tmp_file.write("\t".join(fields) + "\n")

                    close_properly(tmp_file)
                    string = tmp_file.name

        # we will work with string
        if 'w' in self._mode:
            self._mode = 'w'

        return super(FormattedFile, self).__call__(string)


# ---------------------------------------------------------------
# Check all files exist in a glob
# ---------------------------------------------------------------


class globbedFileList(argparse.Action):
    """
    Check the files exist.
    """

    def __init__(self,
                 option_strings,
                 dest,
                 nargs=None,
                 const=None,
                 default=None,
                 type=None,
                 choices=None,
                 required=False,
                 help=None,
                 metavar=None):
        argparse.Action.__init__(self,
                                 option_strings=option_strings,
                                 dest=dest,
                                 nargs=nargs,
                                 const=const,
                                 default=default,
                                 type=type,
                                 choices=choices,
                                 required=required,
                                 help=help,
                                 metavar=metavar,
                                 )

    def __call__(self,
                 parser,
                 namespace,
                 values,
                 option_string=None):

        if "*" in values:
            values = glob.glob(values)
        else:
            values = [values]

        check_file_or_dir_exists(values)

        # values = [ io.TextIOWrapper(gzip.GzipFile(i, 'r')) if ".gz" in i else open(i, "r") for i in values]

        values = [open(i, "r") for i in values]

        # Add the attribute
        setattr(namespace, self.dest, values)

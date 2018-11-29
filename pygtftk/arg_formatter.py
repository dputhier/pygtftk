# -*- coding: utf-8 -*-
"""
Command Line Interface display format.
"""

from __future__ import absolute_import
from __future__ import print_function

import argparse
import glob
import io
import operator
import os
import re
import sys
from builtins import str
from functools import partial

from pybedtools import BedTool

from pygtftk.utils import check_file_or_dir_exists
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
                assert ' '.join(opt_parts) == opt_usage
                assert ' '.join(pos_parts) == pos_usage

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

    def type_func(a_value):
        if val_type == "int":
            try:
                a_value = int(a_value)
            except:
                raise argparse.ArgumentTypeError("Not an int.")
        else:
            try:
                a_value = float(a_value)
            except:
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
# Check a numeric is an int greater than 0
# ---------------------------------------------------------------

def int_greater_than_null(a_value):
    """Check a numeric is an int greater than 0."""
    try:
        a_value = int(a_value)
    except:
        raise argparse.ArgumentTypeError("An integer with min value 1.")
    if a_value <= 0:
        raise argparse.ArgumentTypeError("Minimum value is 1.")
    return a_value


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
        check_file_or_dir_exists(values)
        values = open(values, "r")
        chrom_info_as_dict(values)

        # Add the attribute
        setattr(namespace, self.dest, values)


# ---------------------------------------------------------------
# Returns chromfile as a dict
# ---------------------------------------------------------------


class chromFileAsDict(argparse.Action):
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
        check_file_or_dir_exists(values)
        values = open(values, "r")
        values = chrom_info_as_dict(values)

        # Add the attribute
        setattr(namespace, self.dest, values)


# ---------------------------------------------------------------
# Check file is in bed6 format
# ---------------------------------------------------------------

class bed6(argparse.Action):
    """
    Check the BED file exist and has proper format.
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

        check_file_or_dir_exists(values)

        try:
            file_bo = BedTool(values)
            a = len(file_bo)
        except:
            msg = "Unable to load file: " + values + "."
            message(msg, type="ERROR")
            sys.exit()

        if len(file_bo) == 0:
            msg = "It seems that file " + values + " is empty."
            message(msg, type="ERROR")
            sys.exit()

        if file_bo.file_type != 'bed':
            msg = "File {f} is not a valid bed file."
            msg = msg.format(f=values)
            message(msg, type="ERROR")
            sys.exit()

        names = set()

        for line in file_bo:
            if len(line.fields) != 6:
                message("File -- " + values + " --Need a BED6 file.", type="ERROR")
                sys.exit()

            if line.strand not in ['-', '+', '.']:
                message("File -- " + values + " -- strand is not in proper format", type="ERROR")
            if line.name not in names:
                names.add(line.name)
            else:
                message("File -- " + values + " -- Names (4th columns of the bed file) should be unambiguous.",
                        type="ERROR")
                sys.exit()

        # Add the attribute
        setattr(namespace, self.dest, values)


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


# ---------------------------------------------------------------
# Check file extension
# ---------------------------------------------------------------


class FileWithExtension(argparse.FileType):
    """
    Declare and check file extensions.
    """

    def __init__(self, mode='r', valid_extensions=None, **kwargs):
        super(FileWithExtension, self).__init__(mode, **kwargs)
        self.valid_extensions = valid_extensions

    def __call__(self, string):
        match = False
        if self.valid_extensions:
            if isinstance(self.valid_extensions, str):
                if not string.endswith(self.valid_extensions):
                    if re.search(self.valid_extensions, string):
                        match = True
                else:
                    match = True
            elif isinstance(self.valid_extensions, tuple):
                for exp in self.valid_extensions:
                    if not string.endswith(exp):
                        if re.search(exp, string):
                            match = True
                            break
                    else:
                        match = True
                        break

        if not match:
            message('Not a valid filename extension :' + string, type="WARNING")
            message('Extension expected: ' + str(self.valid_extensions),
                    type="ERROR")
            sys.exit()

        outputdir = os.path.dirname(os.path.abspath(string))

        if not os.path.exists(outputdir):
            if 'w' in self._mode:
                message("Directory not found. Creating.", type="WARNING")
                os.makedirs(outputdir)

        # we will work with string
        if 'w' in self._mode:
            self._mode = 'w'

        return super(FileWithExtension, self).__call__(string)


# ---------------------------------------------------------------
# Pre-defined file types with extension constrains
# ---------------------------------------------------------------


# gtf
gtf_rwb = partial(FileWithExtension, valid_extensions='\.[Gg][Tt][Ff](\.[Gg][Zz])?$')

# gtf file not gz
gtf_rw = partial(FileWithExtension, valid_extensions='\.[Gg][Tt][Ff]$')

# gtf file or txt file
gtf_or_txt_rw = partial(FileWithExtension, valid_extensions=('\.[Gg][Tt][Ff]$',
                                                             '\.[Tt][Xx][Tt]',
                                                             '\.[Cc][Ss][Vv]',
                                                             '\.[Tt][Aa][Bb]',
                                                             '\.[Tt][Ss][Vv]'))

# bed file not gz
bed_rw = partial(FileWithExtension, valid_extensions=('\.[Bb][Ee][Dd]$',
                                                      '\.[Bb][Ee][Dd]3$',
                                                      '\.[Bb][Ee][Dd]6$'))
# txt file
txt_rw = partial(FileWithExtension, valid_extensions=('\.[Tt][Xx][Tt]',
                                                      '\.[Cc][Ss][Vv]',
                                                      '\.[Dd][Ss][Vv]',
                                                      '\.[Tt][Aa][Bb]',
                                                      '\.[Tt][Ss][Vv]'))
# bigwig file
bw_rw = partial(FileWithExtension, valid_extensions=('\.[Bb][Ww]$',
                                                     '\.[Bb][Ii][Gg][Ww][Ii][Gg]$'))

# bigwig file
gtf_or_bed_rwb = partial(FileWithExtension, valid_extensions=('\.[Gg][Tt][Ff](\.[Gg][Zz])?$',
                                                              '\.[Bb][Ee][Dd]$',
                                                              '\.[Bb][Ee][Dd]3$',
                                                              '\.[Bb][Ee][Dd]6$'))

# fasta file
fasta_rw = partial(FileWithExtension, valid_extensions=('\.[Ff][Aa][Ss][Tt][Aa]$',
                                                        '\.[Ff][Nn][Aa]$',
                                                        '\.[Ff][Aa]$',
                                                        '\.[Ff][Aa][Ss]$',
                                                        '\.[Ff][Ff][Nn]$',
                                                        '\.[Ff][Rr][Nn]$'))

# zip file
zip_rw = partial(FileWithExtension, valid_extensions='\.[Zz][Ii][Pp]$')

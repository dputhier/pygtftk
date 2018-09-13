# -*- coding: utf-8 -*-
"""
Command Line Interface display format.
"""

from __future__ import absolute_import
from __future__ import print_function

import argparse
import glob
import os
import re
import sys
from collections import defaultdict

from builtins import object
from builtins import range
from builtins import str
from pybedtools import BedTool

from pygtftk.utils import PY3
from pygtftk.utils import check_file_or_dir_exists
from pygtftk.utils import chrom_info_as_dict
from pygtftk.utils import message

# ---------------------------------------------------------------
# Python2/3  compatibility
# ---------------------------------------------------------------


try:
    basestring
except NameError:
    basestring = str

if PY3:
    from io import IOBase

    file = IOBase


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
                    if isinstance(action.default, file):
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
            format = self._format_actions_usage
            action_usage = format(optionals + positionals, groups)
            usage = ' '.join([s for s in [prog, action_usage] if s])

            # wrap the usage parts if it's too long
            text_width = self._width - self._current_indent
            if len(prefix) + len(usage) > text_width:

                # break usage into wrappable parts
                part_regexp = r'\(.*?\)+|\[.*?\]+|\S+'
                opt_usage = format(optionals, groups)
                pos_usage = format(positionals, groups)
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


class NumericRange(object):
    """To be used in argparse. Ensure float is between start and end.

    :param start: lower range
    :param end: upper range.

    :Example:

    >>> from pygtftk.arg_formatter import NumericRange
    >>> assert 100 == NumericRange(100,200)
    >>> assert 99 != NumericRange(100,200)
    >>> assert 200 == NumericRange(100,200)
    >>> assert 201 != NumericRange(100,200)

    """

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __repr__(self):
        msg = '{0} to {1}'
        msg = msg.format(self.start, self.end)
        return msg

    def __eq__(self, other):
        return self.start <= other <= self.end


class NumericGreaterOrEqual(object):
    """To be used in argparse. Ensure the numeric is greater or equal a define
    value.

    :param start: lower range

    :Example:

    >>> from pygtftk.arg_formatter import NumericGreaterOrEqual
    >>> assert str(NumericGreaterOrEqual(100)) == "100 or more"
    >>> assert 100 == NumericGreaterOrEqual(100)
    >>> assert 200 == NumericGreaterOrEqual(100)
    >>> assert 20 != NumericGreaterOrEqual(100)
    """

    def __init__(self, start):
        self.start = start

    def __repr__(self):
        msg = '{0} or more'
        msg = msg.format(str(self.start))
        return msg

    def int(self):
        return int(self.start)

    def __eq__(self, other):
        return int(other) >= int(self.start)


def int_ge_to_null(a_value):
    """Check a numeric is an int greater or equal to 0."""
    try:
        a_value = int(a_value)
    except:
        raise argparse.ArgumentTypeError("An integer with min value 1.")

    if a_value < 0:
        raise argparse.ArgumentTypeError("Minimum value is 1.")
    return a_value


def int_greater_than_null(a_value):
    """Check a numeric is an int greater than 0."""
    try:
        a_value = int(a_value)
    except:
        raise argparse.ArgumentTypeError("An integer with min value 1.")
    if a_value <= 0:
        raise argparse.ArgumentTypeError("Minimum value is 1.")
    return a_value


def int_greater_than_null_or_None(a_value):
    """Check a numeric is an int greater than 0."""

    if a_value is not None:
        try:
            a_value = int(a_value)
        except:
            raise argparse.ArgumentTypeError("An integer with min value 1.")
        if a_value <= 0:
            raise argparse.ArgumentTypeError("Minimum value is 1.")
    return a_value


def float_greater_than_null(a_value):
    """Check a numeric is a float greater than 0."""
    try:
        a_value = float(a_value)
    except:
        raise argparse.ArgumentTypeError("Should be a float greater than 0.")
    if a_value <= 0:
        raise argparse.ArgumentTypeError("Should be a float greater than 0.")
    return a_value


def float_grt_than_null_and_lwr_than_one(a_value):
    """Check a numeric is a float greater than 0."""
    try:
        a_value = float(a_value)
    except:
        raise argparse.ArgumentTypeError("Should be a float greater than 0.")
    if a_value <= 0:
        raise argparse.ArgumentTypeError("Should be a float greater than 0.")
    if a_value > 1:
        raise argparse.ArgumentTypeError(
            "Should be a float lower or equal to 1.")
    return a_value


class SeparatedList(object):
    """To be used in argparse. Ensure this is a comma separated list of
    the required type and required length.

    :param length: the required length.
    :param type_arg: the required type.
    :param sep: the separator.
    :param check: also check that for two elements l[i] and l[+i] of the splited
    list, 'l[i] op l[+i]' is True. The operator op can be one of ">", "<", "!=",
    ">=" or ">=".

    :Example:

    >>> from pygtftk.arg_formatter import SeparatedList
    >>> assert "1,2" == SeparatedList(length=2, type_arg=int, sep=",", check="<")
    >>> assert "2,1" == SeparatedList(length=2, type_arg=int, sep=",", check=">")
    >>> assert "2,1" == SeparatedList(length=2, type_arg=int, sep=",", check=None)
    >>> assert "1-2" != SeparatedList(length=3, type_arg=int, sep="-")
    """

    def __init__(self,
                 length=2,
                 type_arg=int,
                 sep=",",
                 check=None):
        self.length = length
        self.type_arg = type_arg
        self.sep = sep
        self.check = check

        assert check in [">", "<", "!=", ">=", ">=", None]

    def __repr__(self):
        sep = self.sep
        if sep == " ":
            sep = "space"

        msg = ': a separated list ("{2}") of {1} with {0} element(s)'
        msg = msg.format(str(self.length),
                         re.sub("'.*",
                                "",
                                re.sub("^.*?'",
                                       "",
                                       str(self.type_arg))), sep)
        if self.check is not None:
            msg += ". Constrain: " + "list[i] " + self.check + "list[i + 1]."

        return msg

    def __eq__(self, other):

        assert isinstance(other, basestring)
        other = other.split(self.sep)

        if self.type_arg == int:
            try:
                for i in range(len(other)):
                    other[i] = int(other[i])
            except:
                return False

        if self.check is not None:
            if len(other) > 1:
                for i in range(len(other) - 1):
                    val_1 = other[i]
                    val_2 = other[i + 1]
                    to_test = "val_1 " + self.check + " val_2"

                    if not eval(to_test):
                        return False

        return len(other) == self.length


class NonChromFileError(Exception):
    """
    Raised when we expect a chromosome file and did not get one
    """
    pass


class checkChromFile(argparse.Action):
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


class bedFileList(argparse.Action):
    """
    Check these are bed files.
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

        values = values.split(",")

        for i in values:

            check_file_or_dir_exists(i)

            try:
                file_bo = BedTool(i)
                a = len(file_bo)
            except:
                msg = "Unable to load file: " + i + "."
                message(msg, type="ERROR")
                sys.exit()

            if len(file_bo) == 0:
                msg = "It seems that file " + i.fn + " is empty."
                message(msg, type="ERROR")
                sys.exit()

            if file_bo.file_type != 'bed':
                msg = "File {f} is not a valid bed file."
                msg = msg.format(f=i.fn)
                message(msg, type="ERROR")
                sys.exit()

        values = [open(i, "r") for i in values]

        # Add the attribute
        setattr(namespace, self.dest, values)


class bedFile(argparse.FileType):
    """
    Check the BED file exist and has proper format.
    """

    def __init__(self, **kwargs):

        super(bedFile, self).__init__(**kwargs)

    def __call__(self, string):

        i = string

        check_file_or_dir_exists(i)

        try:
            file_bo = BedTool(i)
            a = len(file_bo)
        except:
            msg = "Unable to load file: " + i.fn + "."
            message(msg, type="ERROR")
            sys.exit()

        if len(file_bo) == 0:
            msg = "It seems that file " + i.fn + " is empty."
            message(msg, type="ERROR")
            sys.exit()

        if file_bo.file_type != 'bed':
            msg = "File {f} is not a valid bed file."
            msg = msg.format(f=i.fn)
            message(msg, type="ERROR")
            sys.exit()

        return super(bedFile, self).__call__(string)


class bedFileWithUnambiguousNames(argparse.FileType):
    """
    Check the BED file exist and has proper format.
    """

    def __init__(self, **kwargs):

        super(bedFileWithUnambiguousNames, self).__init__(**kwargs)

    def __call__(self, string):

        bedfile = string

        check_file_or_dir_exists(bedfile)

        try:
            file_bo = BedTool(bedfile)
            a = len(file_bo)
        except:
            msg = "Unable to load file: " + bedfile + "."
            message(msg, type="ERROR")
            sys.exit()

        if len(file_bo) == 0:
            msg = "It seems that file " + bedfile + " is empty."
            message(msg, type="ERROR")
            sys.exit()

        if file_bo.file_type != 'bed':
            msg = "File {f} is not a valid bed file."
            msg = msg.format(f=bedfile)
            message(msg, type="ERROR")
            sys.exit()

        names = defaultdict(int)

        for line in file_bo:
            if len(line.fields) < 6:
                message("Need a BED6 file.", type="ERROR")
                sys.exit()

            names[line.name] += 1
            if names[line.name] > 1:
                message("File -- " + bedfile + " -- Names (4th columns of the bed file) should be unambiguous.",
                        type="ERROR")
                sys.exit()

        return super(bedFileWithUnambiguousNames, self).__call__(string)


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
        values = glob.glob(values)

        check_file_or_dir_exists(values)

        values = [open(i, "r") for i in values]

        # Add the attribute
        setattr(namespace, self.dest, values)


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
            if isinstance(self.valid_extensions, basestring):
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
            message('Not a valid filename extension', type="WARNING")
            message('Extension expected: ' + str(self.valid_extensions),
                    type="ERROR")
            sys.exit()

        outputdir = os.path.dirname(os.path.abspath(string))

        if not os.path.exists(outputdir):
            if 'w' in self._mode:
                message("Directory not found. Creating.", type="WARNING")
                os.makedirs(outputdir)

        # we will work with string
        if PY3:
            if 'w' in self._mode:
                self._mode = 'w'

        return super(FileWithExtension, self).__call__(string)

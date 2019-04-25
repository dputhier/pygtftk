""" The command manager is intended to store command object
and their associated functions."""

import argparse
import errno
import glob
import io
import logging
import os
import re
import shutil
import subprocess
import sys
import textwrap
from argparse import Action
from importlib.machinery import SourceFileLoader
from sys import platform

import cloudpickle
import yaml

import pygtftk
import pygtftk.cmd_object
import pygtftk.plugins
import pygtftk.settings
import pygtftk.utils
from pygtftk.arg_formatter import ArgFormatter
from pygtftk.arg_formatter import ranged_num
from pygtftk.utils import left_strip_str
from pygtftk.utils import make_tmp_dir
from pygtftk.utils import message
from pygtftk.utils import mkdir_p
from pygtftk.version import __version__


# ---------------------------------------------------------------
# Changing argparse._SubParsersAction to allow subgroups of
# commands.
# https://stackoverflow.com/questions/44292006/display-subgroups-of-commands-with-argparse-subparser
# ---------------------------------------------------------------


class _PseudoGroup(Action):

    def __init__(self, container, title):
        sup = super(argparse._SubParsersAction._PseudoGroup, self)
        sup.__init__(option_strings=[], dest=title)
        self.container = container
        self._choices_actions = []

    def add_parser(self, name, **kwargs):
        # add the parser to the main Action, but move the pseudo action
        # in the group's own list
        parser = self.container.add_parser(name, **kwargs)
        choice_action = self.container._choices_actions.pop()
        self._choices_actions.append(choice_action)
        return parser

    def _get_subactions(self):
        return self._choices_actions

    def add_parser_group(self, title):
        # the formatter can handle recursive subgroups
        grp = argparse._SubParsersAction._PseudoGroup(self, title)
        self._choices_actions.append(grp)
        return grp


argparse._SubParsersAction._PseudoGroup = _PseudoGroup


def add_parser_group(self, title):
    grp = argparse._SubParsersAction._PseudoGroup(self, title)
    self._choices_actions.append(grp)
    return grp


argparse._SubParsersAction.add_parser_group = add_parser_group


# ---------------------------------------------------------------
# An additional action that print Bash completion
# ---------------------------------------------------------------

class BashCompletionAction(argparse._StoreTrueAction):
    """A class to be used by argparser to get bash completion."""

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        super(BashCompletionAction, self).__init__(
            option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        print(pygtftk.settings.get_completion_script())
        sys.exit()


# ---------------------------------------------------------------
# An additional action that list available plugins
# ---------------------------------------------------------------


class ListPlugins(argparse._StoreTrueAction):
    """A class to be used by argparser to list plugins."""

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        super(ListPlugins, self).__init__(
            option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        print("\n".join(list(CmdManager.cmd_obj_list.keys())))
        sys.exit()


# ---------------------------------------------------------------
# An additional action to get all tests
# ---------------------------------------------------------------


class GetTests(argparse._StoreTrueAction):
    """A class to be used by argparser to get all plugin tests
    and write them to a file."""

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        super(GetTests, self).__init__(
            option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):

        cmd_with_no_test = list()

        for cmd in CmdManager.cmd_obj_list:

            test = CmdManager.cmd_obj_list[cmd].test

            if test is not None and re.search("@test", test):
                if platform == "darwin":
                    print(test)
                else:
                    test = test.replace("md5sum-lite", "md5sum")
                    print(test)
            else:
                cmd_with_no_test += [cmd]

        clean = '''

        #clean all
        @test "clean" {
         result=`rm -Rf H3K4me3_cond_* control_list \#H3K4me3_cond_1.bed# *.bw *.bed *.pdf *.png *.genome* *gtf.gz *bed.gz *gtf  *.toto* simpl_* simple_join.txt simple_mat.zip simple_mat simple.2.bw simple.3.bw  ologram_output *_test.bats *simple_mat* simple*.bw`
          [ "$result" = "" ]
        }

        '''

        print(clean)

        for cmd in cmd_with_no_test:
            print("# WARNING: No test found for plugin " + cmd + ".")

        sys.exit()


class GetTestsNoCon(argparse._StoreTrueAction):
    """A class to be used by argparser to get all plugin tests (that do not require internet conn)
    and write them to a file."""

    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        super(GetTestsNoCon, self).__init__(
            option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):

        cmd_with_no_test = list()

        for cmd in CmdManager.cmd_obj_list:

            test = CmdManager.cmd_obj_list[cmd].test

            if cmd not in ['retrieve', 'select_by_go']:
                if test is not None and re.search("@test", test):
                    if platform == "darwin":
                        print(test)
                    else:
                        test = test.replace("md5sum-lite", "md5sum")
                        print(test)
                else:
                    cmd_with_no_test += [cmd]

        clean = '''

        #clean all
        @test "clean" {
         result=`rm -Rf H3K4me3_cond_* control_list \#H3K4me3_cond_1.bed# *.bw *.bed *.pdf *.png *.genome* *gtf.gz *bed.gz *gtf  *.toto* simpl_* simple_join.txt simple_mat.zip simple_mat simple.2.bw simple.3.bw  ologram_output *_test.bats *simple_mat* simple*.bw`
          [ "$result" = "" ]
        }

        '''

        print(clean)

        for cmd in cmd_with_no_test:
            print("# WARNING: No test found for plugin " + cmd + ".")

        sys.exit()


# ---------------------------------------------------------------
# An additional action to add plugins
# ---------------------------------------------------------------


class AddPlugin(argparse.Action):
    """A class to be used by argparser to install plugins."""

    def __init__(self, option_strings, dest, nargs='+', **kwargs):
        super(AddPlugin, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):

        if not os.path.exists(CmdManager.config_dir):
            message("Please run gtftk -h before adding additional plugins",
                    force=True)
            sys.exit()
        values = values.split(" ")[0]
        values = values.split(",")
        plugin_dir_user = os.path.join(CmdManager.config_dir,
                                       "plugins")

        message("Searching for new plugins.", force=True)

        dir_path = make_tmp_dir(prefix="gtftk_AddPlugin")

        dir_path = os.path.join(dir_path, "gtftk")
        cmd = "git clone " + values[0] + " " + dir_path
        p = subprocess.Popen(cmd, shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT)

        while True:
            line = p.stdout.readline()
            if line != '':
                message(line, force=True)
            if line == '' and p.poll() is not None:
                break

        if len(values) == 1:
            f_list = glob.glob(os.path.join(dir_path, "*.py"))
        elif len(values) == 2:
            f_list = glob.glob(os.path.join(dir_path, values[1], "*.py"))
        elif len(values) == 3:
            f_list = glob.glob(os.path.join(dir_path, values[1], values[2]))

        if len(f_list) > 0:

            for f in f_list:
                message("Retrieving plugins %s." %
                        os.path.basename(f),
                        force=True)
                shutil.copy(f, plugin_dir_user)

            open(os.path.join(CmdManager.config_dir, "reload"), "w")
            message("New plugins will be loaded at next startup.",
                    force=True)

        else:
            message("No plugin found.",
                    force=True)

        sys.exit()


# ---------------------------------------------------------------
# An additional action to force plugin reload
# ---------------------------------------------------------------


class UpdatePlugin(argparse._StoreTrueAction):
    """A class to be used by argparser to install plugins."""

    def __init__(self, option_strings, dest, nargs='+', **kwargs):
        super(UpdatePlugin, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        if not os.path.exists(CmdManager.config_dir):
            message("Please run gtftk -h before adding additional plugins",
                    force=True)
            sys.exit(0)

        open(os.path.join(CmdManager.config_dir, "reload"), "w")
        message("Plugins will be updated at next startup.", force=True)
        sys.exit()


# ---------------------------------------------------------------
# An additional action to get information about the environment
# ---------------------------------------------------------------

class GetSysInfo(argparse._StoreTrueAction):
    """A class to be used by argparser to get information about the environment."""

    def __init__(self, option_strings, dest, nargs='+', **kwargs):
        super(GetSysInfo, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        from pandas import __version__ as pandas_ver
        from pybedtools import __version__ as pybedtools_ver
        from pyBigWig import __version__ as bigwig_ver
        from pygtftk import __path__ as pygtftk_path
        import subprocess
        from pygtftk.utils import chomp
        info_sys = []
        info_sys += ['\n- pygtftk version : ' + __version__]
        info_sys += ['- pygtftk path : ' + pygtftk_path[0]]
        info_sys += ['- python version : ' + str(sys.version_info)]
        info_sys += ['- python path : ' + str(sys.prefix)]
        info_sys += ['- pandas version : ' + pandas_ver]
        bedtools_ver = chomp(subprocess.Popen("bedtools --version",
                                              shell=True,
                                              stdout=subprocess.PIPE).stdout.read().decode())
        info_sys += ['- Bedtools version : ' + bedtools_ver]
        info_sys += ['- pybedtools version : ' + pybedtools_ver]
        info_sys += ['- pyBigWig version : ' + bigwig_ver]
        info_sys += ['- uname : ' + str(os.uname())]
        print("\n".join(info_sys))
        sys.exit()


# ---------------------------------------------------------------
# An additional action to get plugin path
# ---------------------------------------------------------------

class GetPluginPath(argparse._StoreTrueAction):
    """A class to be used by argparser to plugin path."""

    def __init__(self, option_strings, dest, nargs='+', **kwargs):
        super(GetPluginPath, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        print(CmdManager.config_dir)
        sys.exit()


# ---------------------------------------------------------------
# The cmdManager class
# ---------------------------------------------------------------


class CmdManager(object):
    """
    This class is intended to store command objects, their associated functions
    and their arguments (that will be exposed through argparse).
    """

    # parsed args

    args = None
    config_dir = None
    config_ffd = None
    config_file = None
    dumped_plugin_path = None
    version_file = None
    reload = False
    # hash should contain the md5 of
    # the path to the gtftk program.
    # the bin/gtftk program is suppose to
    # fill this attribute
    hash = None

    # -------------------------------------------------------------------------
    # The main parser
    # -------------------------------------------------------------------------

    prg_desc = """
  A toolbox to handle GTF files.

  Example:

  gtftk get_example -f chromInfo -o simple.chromInfo ; 
  gtftk get_example  | gtftk feature_size -t mature_rna | gtftk nb_exons |\\
  gtftk intron_sizes | gtftk exon_sizes | gtftk convergent -u 24 -d 24  -c simple.chromInfo | \\
  gtftk divergent -u 101 -d 10  -c simple.chromInfo  | \\
  gtftk overlapping -u 0 -d 0 -t transcript -c simple.chromInfo -a |  \\
  gtftk select_by_key -k feature -v transcript |   gtftk tabulate -k "*" -b -x
  

  Type 'gtftk sub-command -h' for more information.

    """

    parser = argparse.ArgumentParser(
        formatter_class=ArgFormatter,
        description=prg_desc,
        epilog="------------------------\n"
    )

    parser._optionals.title = "Main command arguments"

    parser.add_argument('-b', '--bash-comp',
                        nargs=0,
                        help="Get a script to activate bash completion.",
                        action=BashCompletionAction)

    parser.add_argument('-p', '--plugin-tests',
                        nargs=0,
                        help="Display bats tests for all plugin.",
                        action=GetTests)

    parser.add_argument('-u', '--plugin-tests-no-conn',
                        nargs=0,
                        help="Display bats tests for plugins not relying on server conn.",
                        action=GetTestsNoCon)

    parser.add_argument('-s', '--system-info',
                        nargs=0,
                        help="Display some info about the system.",
                        action=GetSysInfo)

    parser.add_argument('-d', '--plugin-path',
                        nargs=0,
                        help="Print plugin path",
                        action=GetPluginPath)

    parser.add_argument('-v', '--version',
                        action='version',
                        version='%(prog)s v{0}'.format(__version__))

    parser.add_argument('-l', '--list-plugins',
                        nargs=0,
                        help="Get the list of plugins.",
                        action=ListPlugins)

    # -------------------------------------------------------------------------
    # The sub parser
    # -------------------------------------------------------------------------

    # Declare a subparser
    sub_parsers = parser.add_subparsers(
        title='Available sub-commands/plugins',
        dest='command',
        metavar='')

    # -------------------------------------------------------------------------
    # Declare subparser groups
    # -------------------------------------------------------------------------

    # Declare subparser groups

    # geno_info, edition, selection, conversion, annotation, info, coverage, sequence
    grp_editing = sub_parsers.add_parser_group('\n------- Editing --------\n')
    grp_info = sub_parsers.add_parser_group('\n----- Information ------\n')
    grp_select = sub_parsers.add_parser_group('\n------ Selection -------\n')
    grp_convert = sub_parsers.add_parser_group('\n------ Conversion ------\n')
    grp_annot = sub_parsers.add_parser_group('\n------ Annotation ------\n')
    grp_seq = sub_parsers.add_parser_group('\n------- Sequence -------\n')
    grp_coord = sub_parsers.add_parser_group('\n----- Coordinates ------\n')
    grp_cov = sub_parsers.add_parser_group('\n------- Coverage -------\n')
    grp_misc = sub_parsers.add_parser_group('\n----- Miscellaneous ----\n')

    # -----------------------------------------------------------------------
    # A dict of cmdObjects (plugins)
    # -----------------------------------------------------------------------

    # This class attributes stores the instances of CmdObject
    cmd_obj_list = dict()

    # -------------------------------------------------------------------------
    # Methods
    # -------------------------------------------------------------------------

    @classmethod
    def check_config_file(cls):

        # ----------------------------------------------------------------------
        # Config directory and config files
        # ----------------------------------------------------------------------

        CmdManager.config_dir = os.path.join(os.path.expanduser("~"),
                                             ".gtftk",
                                             CmdManager.hash)

        CmdManager.config_file = os.path.join(CmdManager.config_dir,
                                              "gtftk.cnf")

        CmdManager.dumped_plugin_path = os.path.join(CmdManager.config_dir,
                                                     "plugin.pick")

        CmdManager.version_file = os.path.join(CmdManager.config_dir,
                                               "version.py")

        if os.path.exists(os.path.join(CmdManager.config_dir, "reload")):

            CmdManager.reload = True
            os.unlink(os.path.join(CmdManager.config_dir, "reload"))
        else:
            CmdManager.reload = False

        # ----------------------------------------------------------------------
        # Load config
        # ----------------------------------------------------------------------

        # Load the config file  (~/.pygtftk/pygtftk.conf)
        # or create it.
        # This file contains the path to the directories containing plugins.
        # The plugin_path is set to
        # ~/.pygtftk/plugins by default

        if not os.path.exists(CmdManager.config_dir):
            try:
                os.makedirs(CmdManager.config_dir)
            except OSError as exc:
                if exc.errno != errno.EEXIST:
                    raise
                pass

        plug_dir_default = os.path.join(CmdManager.config_dir, "plugins")

        if not os.path.exists(plug_dir_default):
            try:
                os.makedirs(plug_dir_default)
            except OSError as exc:
                if exc.errno != errno.EEXIST:
                    raise
                pass

        if not os.path.exists(CmdManager.config_file):
            with open(CmdManager.config_file, 'w') as a_file:
                a_file.write("---\n")
                out_dict = {'plugin_path': os.path.join(CmdManager.config_dir,
                                                        'plugins')}
                a_file.write(yaml.dump(out_dict, default_flow_style=False))

        # ----------------------------------------------------------------------
        # Check version
        # ----------------------------------------------------------------------

        if os.path.exists(CmdManager.version_file):

            cur_version = None

            version_fh = open(CmdManager.version_file, "r")

            for i in version_fh:
                if "__version__" in i:
                    cur_version = i.split("=")[1]
                    cur_version = re.sub("[\'\" \n\r]", "", cur_version)

            if cur_version != __version__:
                cls._create_version_file()
                CmdManager.reload = True
        else:
            cls._create_version_file()

    @classmethod
    def _create_version_file(cls):
        version_file_installed = os.path.join(pygtftk.__path__[0], "version.py")
        shutil.copy(version_file_installed, CmdManager.config_dir)

    def __init__(self):
        """The constructor."""

        self.check_config_file()

    @classmethod
    def add_command(cls, cmd):
        """ Add new argument parser for a command. Note that
        verbosity, keep-temp and help arguments are added by default."""

        # ----------------------------------------------------------------------
        # Command help display
        # ----------------------------------------------------------------------

        cmd.desc = "  Description: \n     *" + textwrap.fill(
            textwrap.dedent(
                left_strip_str(
                    cmd.desc)).strip(),
            100, initial_indent='  ', subsequent_indent='     ')

        if cmd.notes is not None:
            cmd.notes = cmd.notes.lstrip("\n")
            cmd.notes = cmd.notes.strip()
            cmd.notes = cmd.notes.strip("--")
            cmd.notes = cmd.notes.strip(" ")
            tokens = cmd.notes.split("--")
            cmd.desc += "\n\n" + textwrap.fill(
                textwrap.dedent("Notes:").strip(),
                100, initial_indent='  ', subsequent_indent='     ')

            for i in range(len(tokens)):
                cmd.desc += "\n" + textwrap.fill(
                    textwrap.dedent(
                        " * " +
                        left_strip_str(
                            tokens[i].strip())).strip(),
                    100, initial_indent='     ', subsequent_indent='     ')

        if cmd.references is not None:
            cmd.references = cmd.references.lstrip("\n")
            cmd.references = cmd.references.strip()
            cmd.references = cmd.references.strip("--")
            cmd.references = cmd.references.strip(" ")
            tokens = cmd.references.split("--")
            cmd.desc += "\n\n" + textwrap.fill(
                textwrap.dedent("References:").strip(),
                100, initial_indent='  ', subsequent_indent='     ')

            for i in range(len(tokens)):
                cmd.desc += "\n" + textwrap.fill(
                    textwrap.dedent(
                        " * " +
                        left_strip_str(
                            tokens[i].strip())).strip(),
                    100, initial_indent='     ', subsequent_indent='     ')

        cmd.desc = cmd.desc + "\n\n" + textwrap.fill(
            textwrap.dedent(
                "  Version: " +
                left_strip_str(
                    cmd.updated)).strip(),
            100, initial_indent='  ', subsequent_indent='     ')

        cmd.desc = re.sub("\-\\\-", "--", cmd.desc)

        # ----------------------------------------------------------------------
        # Define command-wise args
        # ----------------------------------------------------------------------

        if cmd.lang == "Python":
            group = cmd.parser.add_argument_group(
                'Command-wise optional arguments')

        # help is a default argument for any command
        group.add_argument("-h",
                           "--help",
                           action="help",
                           help="Show this help message and exit.")

        # verbose is a default argument of any command
        group.add_argument("-V",
                           "--verbosity",
                           default=0,
                           type=ranged_num(lowest=0,
                                           highest=3,
                                           val_type="int",
                                           linc=True,
                                           hinc=True),
                           help="Set output verbosity ([0-3]).",
                           nargs='?',
                           required=False)

        # verbose is a default argument of any command
        group.add_argument("-D",
                           "--no-date",
                           action="store_true",
                           help="Do not add date to output file names.")

        # verbose is a default argument of any command
        group.add_argument("-C",
                           "--add-chr",
                           action="store_true",
                           help="Add 'chr' to chromosome names before printing output.")

        # keep-temp-file is a default argument of any command
        group.add_argument("-K",
                           "--tmp-dir",
                           type=str,
                           default=None,
                           help="Keep all temporary files into this folder.",
                           required=False)

        # keep-temp-file is a default argument of any command
        group.add_argument("-A",
                           "--keep-all",
                           action="store_true",
                           help="Try to keep all temporary files even if process does not terminate normally.",
                           required=False)

        # logger-file can be used to store the requested command
        # arguments into a file.
        group.add_argument("-L",
                           "--logger-file",
                           type=str,
                           help='Stores the arguments passed to the command into a file.',
                           required=False)

        # logger-file can be used to store the requested command
        # arguments into a file.
        group.add_argument("-W",
                           "--write-message-to-file",
                           type=argparse.FileType('w'),
                           default=None,
                           help='Store all message into a file.',
                           required=False)

        # Add the command to the list of known command
        cls.cmd_obj_list[cmd.name] = cmd

        # Update the global argument parser

        arg_dict = {'name': cmd.name,
                    'formatter_class': ArgFormatter,
                    'parents': [cmd.parser],
                    'help': cmd.message,
                    'add_help': False,
                    'description': cmd.desc}
        if cmd.group == 'editing':
            cls.grp_editing.add_parser(**arg_dict)

        elif cmd.group == 'information':
            cls.grp_info.add_parser(**arg_dict)

        elif cmd.group == 'selection':
            cls.grp_select.add_parser(**arg_dict)

        elif cmd.group == 'conversion':
            cls.grp_convert.add_parser(**arg_dict)

        elif cmd.group == 'coordinates':
            cls.grp_coord.add_parser(**arg_dict)

        elif cmd.group == 'annotation':
            cls.grp_annot.add_parser(**arg_dict)

        elif cmd.group == 'sequences':
            cls.grp_seq.add_parser(**arg_dict)

        elif cmd.group == 'coverage':
            cls.grp_cov.add_parser(**arg_dict)

        elif cmd.group == 'miscellaneous':
            cls.grp_misc.add_parser(**arg_dict)

        else:
            raise ValueError("Unknow group for command : %s" % cmd.name)

    @staticmethod
    def _find_plugins():

        message("Searching plugins", force=True)
        config_file = CmdManager.config_file

        # User plugins
        plugin_dir_user = yaml.load(open(config_file, "r"), Loader=yaml.FullLoader)["plugin_path"]
        sys.path.append(plugin_dir_user)
        plugins = sorted(os.listdir(plugin_dir_user))
        plugins_user = [os.path.join(plugin_dir_user, x) for x in plugins]

        # System wide plugins (those declared in the plugins directory of
        # pygtftk source)

        plugin_dir_base = os.path.join(pygtftk.__path__[0], "plugins")
        sys.path.append(plugin_dir_base)
        plugins = sorted(os.listdir(plugin_dir_base))
        plugins_system = [os.path.join(plugin_dir_base, x) for x in plugins]

        plugins = plugins_user + plugins_system

        for plug in plugins:
            if plug.endswith(".py") and plug != "__init__.py":

                # Loading the plugin should force code to create
                # a cmdObject that will be added to the CmdManager
                # gtftk.plugins.tss_dist
                module_name = re.sub("\.py$", "", plug)
                module_name = re.sub("/", ".", module_name)
                module_name = re.sub(".*pygtftk", "pygtftk", module_name)

                try:

                    SourceFileLoader(module_name, plug).load_module()

                except Exception as e:
                    message("Failed to load plugin :" + plug, type="WARNING")
                    print(e)

            elif plug.endswith(".R"):
                pass
                # declare_r_cmd(plugin_path, plug)

        CmdManager.reload = False

        if os.path.exists(os.path.join(CmdManager.config_dir, "reload")):
            os.remove(os.path.join(CmdManager.config_dir, "reload"))

    def dump_plugins(self):
        """Save the plugins into a pickle object."""

        message("Dumping plugins", force=True)

        f_handler = open(CmdManager.dumped_plugin_path, "wb")

        pick = cloudpickle.CloudPickler(f_handler)
        pick.dump((self.cmd_obj_list, self.parser))
        f_handler.close()
        # self.load_plugins()

    def load_plugins(self):
        """Load the plugins."""

        if CmdManager.reload:

            shutil.rmtree(CmdManager.config_dir, ignore_errors=True)
            self.check_config_file()
            self._find_plugins()
            self.dump_plugins()
        else:
            if not os.path.exists(CmdManager.dumped_plugin_path):
                shutil.rmtree(CmdManager.config_dir, ignore_errors=True)
                self.check_config_file()
                self._find_plugins()
                self.dump_plugins()

        self._load_dumped_plugins()

    @staticmethod
    def _load_dumped_plugins():

        f_handler = open(CmdManager.dumped_plugin_path, "rb")
        CmdManager.cmd_obj_list, CmdManager.parser = cloudpickle.load(f_handler)
        f_handler.close()

        for cur_cmd in sorted(CmdManager.cmd_obj_list):

            # fix some issues related to dumping of the parser
            for cur_arg in CmdManager.cmd_obj_list[cur_cmd].parser._option_string_actions:

                obj = CmdManager.cmd_obj_list[
                    cur_cmd].parser._option_string_actions[cur_arg]

                if obj.default == '==stdin==':
                    CmdManager.cmd_obj_list[cur_cmd].parser._option_string_actions[
                        cur_arg].default = sys.stdin
                if obj.default == '==SUPPRESS==':
                    CmdManager.cmd_obj_list[cur_cmd].parser._option_string_actions[
                        cur_arg].default = argparse.SUPPRESS

        CmdManager.parser._option_string_actions['-h'].default = argparse.SUPPRESS
        CmdManager.parser._option_string_actions['--help'].default = argparse.SUPPRESS
        CmdManager.parser._option_string_actions['-v'].default = argparse.SUPPRESS
        CmdManager.parser._option_string_actions['--version'].default = argparse.SUPPRESS

        f_handler.close()

    @classmethod
    def parse_cmd_args(cls):
        """ Parse arguments of all declared commands."""

        CmdManager.args = cls.parser.parse_args(None)
        args = CmdManager.args
        cmd_name = args.command

        if cmd_name is None:
            message("Please provide a subcommand or argument (e.g. -h)", type="WARNING", force=True)
            CmdManager.parser.print_help()
            exit(0)

        lang = cls.cmd_obj_list[cmd_name].lang

        if lang == 'Python':
            if args.tmp_dir is not None:

                if not os.path.exists(args.tmp_dir):
                    msg = "Creating directory {d}."
                    message(msg.format(d=args.tmp_dir), type="INFO")
                    mkdir_p(args.tmp_dir)
                if not os.path.isdir(args.tmp_dir):
                    msg = "{d} is not a directory."
                    message(msg.format(d=args.tmp_dir), type="ERROR")

                pygtftk.utils.TMP_DIR = args.tmp_dir

        return args

    @classmethod
    def run(cls, args):
        """Run the selected command.
        """
        #  Retrieve the ad hoc function
        args = dict(args.__dict__)

        cmd_ob = cls.cmd_obj_list[args['command']]

        # Add a logger to the command object
        cls.cmd_obj_list[args['command']].logger = logging.getLogger(__name__)

        # .pyc -> .py
        fun_path = cmd_ob.fun.rstrip("c")
        tmp_module = SourceFileLoader('tmp_module', fun_path).load_module()
        fun = getattr(tmp_module, args['command'])

        # Save args to log file

        if cmd_ob.lang == "Python":

            # Add 'chr' to the chromosome names
            if args['add_chr']:
                pygtftk.utils.ADD_CHR = 1

            if args['write_message_to_file'] is not None:
                pygtftk.utils.MESSAGE_FILE = args['write_message_to_file']

            if args['logger_file'] is not None:
                if os.path.isdir(args['logger_file']):
                    message("ERROR --logger-file is a directory.",
                            type="ERROR")

                if not os.path.exists(args['logger_file']):
                    logger_file_h = open(args['logger_file'], 'w+')
                    logger_file_h.close()

                log_format = "-->> %(asctime)s - %(name)s - " + \
                             "%(levelname)s - %(message)s"
                datefmt = '%Y-%m-%d %H:%M:%S'

                logging.basicConfig(filename=args['logger_file'],
                                    level=logging.INFO,
                                    format=log_format,
                                    datefmt=datefmt)

                cmd_ob.logger.info("Command: " + " ".join(sys.argv))
                cmd_ob.logger.info("Argument: " + 'command=' + args['command'])

                for key, value in list(args.items()):
                    if isinstance(value, io.IOBase):
                        value = value.name
                    else:
                        value = str(value)

                    cmd_ob.logger.info("Argument: " + key + "=" + value)

            # Set the level of verbosity
            # Can be None if -V is used without value
            # (nargs)
            if args['verbosity'] is None:
                pygtftk.utils.VERBOSITY = 1
            else:
                pygtftk.utils.VERBOSITY = int(args['verbosity'])

            # Set whether date should be added to
            # output file
            if args['no_date']:
                pygtftk.utils.ADD_DATE = False

            arg_list_to_del = ['bash_comp', 'add_chr', 'version', 'help',
                               'plugin_tests', 'list_plugins', 'plugin_tests_no_conn',
                               'r_libs', 'add_plugin', 'update_plugins', 'system_info',
                               'plugin_path', 'no_date', 'keep_all', 'logger_file',
                               'tmp_dir', 'verbosity', 'command', 'write_message_to_file']
            for cur_arg in arg_list_to_del:
                try:
                    del args[cur_arg]
                except KeyError:
                    pass

            # Run the command
            fun(**args)

        else:
            raise ValueError("Unknow language.")

""" A container for a command."""

import argparse
import re
import sys

import pygtftk
import pygtftk.cmd_manager


class CmdObject(object):
    """ A simple container for a command."""

    def __init__(self,
                 name=None,
                 message=None,
                 parser=None,
                 fun=None,
                 desc=None,
                 lang="Python",
                 notes=None,
                 references=None,
                 updated="",
                 group=None,
                 test=None,
                 rlib=None):

        self.name = name
        self.notes = notes
        self.references = references
        self.group = group
        self.updated = updated
        self.message = message
        for i in parser._option_string_actions:
            if parser._option_string_actions[i].default is sys.stdin:
                parser._option_string_actions[i].default = '==stdin=='

        for i in range(len(parser.__dict__['_actions'])):
            if isinstance(
                    parser.__dict__['_actions'][i], argparse._HelpAction):
                parser.__dict__['_actions'].pop(i)
                break

        try:
            del parser.__dict__['_option_string_actions']['--help']
            del parser.__dict__['_option_string_actions']['-h']
        except:
            pass

        self.parser = parser
        self.fun = fun
        self.desc = desc
        self.logger = None
        self.lang = lang
        self.test = test
        self.rlib = rlib

        if re.search("@test", self.test):
            pygtftk.cmd_manager.CmdManager.add_command(self)
        else:
            pygtftk.utils.message(
                "%s command has no test and won't be installed." % name,
                type="WARNING")

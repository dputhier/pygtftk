"""
This module is intended to retrieve argument parser declared in *.R files.
"""


def declare_r_cmd(plugin_dir, plugin_fn):
    """Take a *.R file as input.
    This file must contain a parser construct using a call to library(argparse).
     The declare_R_cmd function will execute the code that won't be run
    (options(run.main=FALSE)). It will thus only retrieve the code the python
    code that is produced by library(argparse). Indeed, R argparse version is
    just an interface to python. The code will serve to declare a
    novel command.

    :plugin_dir: the directory containing the plugin.
    :plugin_fn: the name of the plugin file (*.R).

    """

    import os
    from pygtftk.cmd_object import CmdObject
    # import rpy2.robjects as robjects

    file_name = os.path.join(plugin_dir, plugin_fn)

    # -----------------------------------------
    # Interpret R code using rpy2
    # and retrieve the python code produced
    # by argparse (R side).
    # -----------------------------------------

    r_instruc = """
                options(run.main=FALSE)
                source('{fn}')
                """
    r_instruc = r_instruc.format(fn=file_name)

    # Ask R to source R code

    robjects.r['options'](warn=-1)
    robjects.r(r_instruc)

    # Retrieve the parser from R

    parser = robjects.r['parser']

    # get the corresponding python code

    imported_parser = "\n".join(list(parser['python_code']))

    parser = compile(imported_parser, '<string>', 'exec')
    exec(parser)
    parser.add_help = False

    # Declare a novel command to the cmdManager
    CmdObject(name=os.path.splitext(plugin_fn)[0],
              message=parser.__dict__["description"],
              parser=parser,
              fun=file_name,
              desc=parser.__dict__["description"],
              lang="R")


def _main():
    declare_r_cmd("../plugins/", "read_file.R")


if __name__ == "__main__":
    _main()

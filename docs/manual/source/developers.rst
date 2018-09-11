Writing your own commands
=============================


Maybe you will be interested in helping us in the future by sharing your own commands. We provide an example below of the development of a very basic command that prints a GTF file.

Just do the following things to declare a new command named **'print_gtf'**:

- Write a script using the template provided below. This script can be run as a stand-alone script or as a gtftk plugin.
- Add this script to the gtftk/plugins or ~/.gtftk/plugins/ directories.
- Call *gtftk -u* to update the plugin list. A new command should be available upon *gtftk -h*.
- You can propose new commands by sending us pull requests.


.. code-block:: python

   #!/usr/bin/env python
   """
   Description: This module is intented to print a GTF.
   Developer: D. Puthier.
   Last modifications: 09 Mar 2016.
   Version: {v}
   """
   
   import sys
   import argparse
   from pygtftk.cmd_object import CmdObject
   from pygtftk.gtf_interface import GTF
   from pygtftk.arg_formatter import FileWithExtension

   #-------------------------------------------------------------------------
   # NOTES: Just place this file into ~/.gtftk/plugins
   # and ask for: 
   #    gtftk --update-plugins
   # gtftk will load the plugin next time it is called.
   #-------------------------------------------------------------------------
   
      
   
   #-------------------------------------------------------------------------
   # Message can be of type "INFO" (default), "WARNING" or "ERROR"
   # When using message it will adapt to the user-declared level of
   # verbosity
   #-------------------------------------------------------------------------
   
   from pygtftk.utils import message

   #-------------------------------------------------------------------------
   # Any temporary file created using the make_tmp_file function
   # will be deleted upon exit or may be kept into a target directory (see -K
   # command-wise argument)
   # This command should thus be used to create any temporary file.
   #-------------------------------------------------------------------------
   
   from pygtftk.utils import make_tmp_file
    
   #-------------------------------------------------------------------------
   # Command information
   #-------------------------------------------------------------------------
   
   __doc__ = """
   This is the doc about the command that will appear when gtftk my_command -h 
   is called...
   """
   
   
   __notes__ = """
   -- A note that will appear when 'gtftk my_command -h' will be called.
   -- Another note. If you want to refer to long form arguments use '\'. e.g -\-distance.
   """
   
   
   #-------------------------------------------------------------------------
   # First define the function/command arguments.
   # Note that the syntax is the same that would be used for a regular program
   # implementing an argument parser.
   # Make use as possible of argparse.FileType and more complexes types as
   # found in gtftk.arg_formatter.
   #-------------------------------------------------------------------------
   
   def make_parser():
      parser = argparse.ArgumentParser(add_help=True)
   
      parser_grp = parser.add_argument_group('Arguments')
   
      parser_grp.add_argument('-i', '--inputfile',
                              help="Path to the GTF file. Default to STDIN",
                              default=sys.stdin,
                              metavar="GTF",
                              type=FileWithExtension('r',
                                                        valid_extensions=('\.[Gg][Tt][Ff](\.[Gg][Zz])?$')))
         
            
      parser_grp.add_argument('-o', '--outputfile',
                              help="Output file.",
                              default=sys.stdout,
                              metavar="GTF",
                              type=FileWithExtension('w',
                                                     valid_extensions=('\.[Gg][Tt][Ff]$')))
         
       return parser
   
   #-------------------------------------------------------------------------
   # Now we declare a main function, as would be done
   # for a regular program
   #-------------------------------------------------------------------------
   
   
   # NB: The verbosity, tmp_dir=None and logger_file are mandatory arguments
   
   def print_gtf(inputfile=None,
                 outputfile=None,
                 tmp_dir=None,
                 logger_file=None,
                 verbosity=0):
       """This function will only print a GTF..."""
   
       message("Reading GTF")
       gtf = GTF(inputfile)
       gtf.write(outputfile)
   
   #-------------------------------------------------------------------------
   # Now we check if the python interpreter is running this module
   # as the main program or whether it is called by the plugin manager.
   #-------------------------------------------------------------------------
   
   def main():
       """The main function."""
       args = make_parser().parse_args()
       args = dict(args.__dict__)
       print_gtf(**args)
       
   if __name__ == '__main__':
       main()   
   else:
   
       # Just declare a new command object
       # That will call the command manager.
       # With the user-passed arguments.
       # Available groups are: editing, information, selection, conversion, 
       # coordinates, annotation, sequences, coverage,
       # and miscellaneous.
   
       cmd = CmdObject(name="print_gtf",
                       message="Print a GTF",
                       parser=make_parser(),
                       fun=os.path.abspath(__file__),
                       group="miscellaneous",
                       desc=__doc__,
                       notes=__notes__)




    
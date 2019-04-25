#!/usr/bin/env python
"""
 Delete transcripts containing an intron whose size is below s.
"""

import argparse
import os
import sys
from collections import OrderedDict
from collections import defaultdict

from pygtftk import arg_formatter
from pygtftk.cmd_object import CmdObject
from pygtftk.gtf_interface import GTF
from pygtftk.utils import close_properly
from pygtftk.utils import message

__updated__ = "2018-01-20"

__notes__ = """
 -- By default delete transcripts containing an intron (or a sum intronic nucleotides, see -m) whose size is below s.
 -- If -\-invert-match is selected delete transcripts containing an intron (or a sum intronic nucleotides, see -m) whose size is greater or equal to s.
 -- Mono-exonic transcripts are not tested for intron size. They can be kept or deleted based on -d.
"""


def make_parser():
    """The program parser."""
    parser = argparse.ArgumentParser(add_help=True)

    parser_grp = parser.add_argument_group('Arguments')

    parser_grp.add_argument('-i', '--inputfile',
                            help="Path to the GTF file. Default to STDIN",
                            default=sys.stdin,
                            metavar="GTF",
                            type=arg_formatter.FormattedFile(mode='r', file_ext=('gtf', 'gtf.gz')))

    parser_grp.add_argument('-o', '--outputfile',
                            help="Output file.",
                            default=sys.stdout,
                            metavar="GTF",
                            type=arg_formatter.FormattedFile(mode='w', file_ext=('gtf')))

    parser_grp.add_argument('-s', '--intron-size',
                            help="The minimum intron size.",
                            type=int,
                            default=100,
                            required=True)

    parser_grp.add_argument('-m', '--merged',
                            help="If selected, the sum of all intron lengths for a gene should be higher than s.",
                            action="store_true")

    parser_grp.add_argument('-d', '--delete-monoexonic',
                            help="Delete mono-exonic transcripts.",
                            action="store_true")

    parser_grp.add_argument('-a', '--add-intron-size',
                            help="Add a new key containing intron_size (comma-separated in order of apppearance) or the sum of intron size (see -m).",
                            action="store_true")

    parser_grp.add_argument('-v', '--invert-match',
                            help='Keep genes with an intron whose size is above s and delete others.',
                            action="store_true",
                            required=False)
    return parser


def select_by_intron_size(
        inputfile=None,
        outputfile=None,
        intron_size=0,
        merged=False,
        invert_match=False,
        delete_monoexonic=False,
        add_intron_size=False):
    """
    Select genes which contain an intron of size at least s or whose sum of intron size is at least s
    """

    message("Searching for intronic regions.")

    gtf = GTF(inputfile, check_ensembl_format=False)

    introns_bo = gtf.get_introns(by_transcript=True,
                                 name=["transcript_id"],
                                 intron_nb_in_name=False).sort()

    # Get the list of transcripts
    all_tx_ids = gtf.get_tx_ids(nr=True)

    # The list of transcripts
    # to be deleted
    to_delete = OrderedDict()

    if merged:
        # Create a dict that will contain the sum of introns for
        # each transcript
        intron_sum_dict = OrderedDict.fromkeys(all_tx_ids, 0)

        for i in introns_bo:
            size = i.end - i.start
            tx_id = i.name
            intron_sum_dict[tx_id] += size

        for tx_id, sum_intron in list(intron_sum_dict.items()):

            if sum_intron != 0:
                if not invert_match:
                    if sum_intron < intron_size:
                        to_delete[tx_id] = 1

                else:
                    if sum_intron >= intron_size:
                        to_delete[tx_id] = 1
            else:
                if delete_monoexonic:
                    to_delete[tx_id] = 1

        if add_intron_size:
            gtf = gtf.add_attr_from_dict(feat="transcript",
                                         key="transcript_id",
                                         a_dict=intron_sum_dict,
                                         new_key="intron_size_sum")

    else:

        # Create a dict that will contain a list introns size
        # for each transcript

        intron_size_dict = defaultdict(list)

        for tx_id in all_tx_ids:
            intron_size_dict[tx_id] = []

        for i in introns_bo:
            size = i.end - i.start
            tx_id = i.name

            intron_size_dict[tx_id] += [size]

        for tx_id, list_size in list(intron_size_dict.items()):
            if not list_size:
                intron_size_dict[tx_id] = [0]
                if delete_monoexonic:
                    to_delete[tx_id] = 1
                continue

            for size in intron_size_dict[tx_id]:

                if not invert_match:
                    if size < intron_size:
                        to_delete[tx_id] = 1

                else:
                    if size >= intron_size:
                        to_delete[tx_id] = 1

        if add_intron_size:

            for tx_id, list_size in list(intron_size_dict.items()):
                list_size = [str(x) for x in list_size]
                intron_size_dict[tx_id] = ",".join(list_size)

            gtf = gtf.add_attr_from_dict(feat="transcript",
                                         key="transcript_id",
                                         a_dict=intron_size_dict,
                                         new_key="intron_size")

    all_tx_ids = gtf.get_tx_ids(nr=True)
    all_tx_ids = [x for x in all_tx_ids if x not in to_delete]
    msg_list = ",".join(list(to_delete.keys()))
    nb_char = min([len(msg_list), 40])
    msg_list = msg_list[0:nb_char]
    message("Deleting: " + msg_list + "...")

    gtf = gtf.select_by_key("transcript_id",
                            ",".join(all_tx_ids))

    gtf.write(outputfile, gc_off=True)

    close_properly(outputfile, inputfile)


def main():
    """The main program."""
    myparser = make_parser()
    args = myparser.parse_args()
    args = dict(args.__dict__)
    select_by_intron_size(**args)


if __name__ == '__main__':
    main()

else:

    test = '''
    #Keep all transcript with intron_size < 1. Monoexonic are kept by default.
    @test "select_by_intron_size_1" {
     result=`gtftk get_example -d simple_04 |  gtftk select_by_intron_size -s 1 -a  -v | gtftk select_by_key -t | wc -l `
      [ "$result" -eq 8 ]
    }

    #Keep all transcript with intron_size < 5. Monoexonic are kept by default.
    @test "select_by_intron_size_2" {
     result=`gtftk get_example -d simple_04 |  gtftk select_by_intron_size -s 5 -a  -v | gtftk select_by_key -t | wc -l `
      [ "$result" -eq 13 ]
    }

    #Keep all transcript with intron_size < 5. Monoexonic are deleted.
    @test "select_by_intron_size_3" {
     result=`gtftk get_example -d simple_04 |  gtftk select_by_intron_size -s 5 -a -d  -v | gtftk select_by_key -t | wc -l `
      [ "$result" -eq 5 ]
    }
        
    #Keep all transcript with intron_size >= 2. Monoexonic are deleted.
    @test "select_by_intron_size_4" {
     result=`gtftk get_example -d simple_04 |  gtftk select_by_intron_size -s 2 -a -d  | gtftk select_by_key -t | wc -l `
      [ "$result" -eq 8 ]
    }

    #Keep all transcript with intron_size >= 2. Monoexonic are kept.
    @test "select_by_intron_size_5" {
     result=`gtftk get_example -d simple_04 |  gtftk select_by_intron_size -s 2 -a   | gtftk select_by_key -t | wc -l `
      [ "$result" -eq 16 ]
    }

    #Keep all transcript with merged intron_size < 1. Monoexonic are kept by default.
    @test "select_by_intron_size_6" {
     result=`gtftk get_example -d simple_04 |  gtftk select_by_intron_size -s 1 -a  -v -m | gtftk select_by_key -t | wc -l `
      [ "$result" -eq 8 ]
    }

    #Keep all transcript with merged intron_size < 1. Monoexonic are deleted. Nothing left
    @test "select_by_intron_size_7" {
     result=`gtftk get_example -d simple_04 |  gtftk select_by_intron_size -s 1 -a  -v -m -d | gtftk select_by_key -t | wc -l `
      [ "$result" -eq 0 ]
    }


    #Keep all transcript with merged intron_size < 5. Monoexonic are deleted.
    @test "select_by_intron_size_8" {
     result=`gtftk get_example -d simple_04 |  gtftk select_by_intron_size -s 5 -m -a -d  -v | gtftk select_by_key -t | wc -l `
      [ "$result" -eq 4 ]
    }

    #Keep all transcript with merged intron_size < 5. Monoexonic are kept.
    @test "select_by_intron_size_3" {
     result=`gtftk get_example -d simple_04 |  gtftk select_by_intron_size -s 5 -m -a  -v | gtftk select_by_key -t | wc -l `
      [ "$result" -eq 12 ]
    }
            
    #Keep all transcript with merged intron_size >= 2. Monoexonic are deleted.
    @test "select_by_intron_size_9" {
     result=`gtftk get_example -d simple_04 |  gtftk select_by_intron_size -s 2 -a -d -m | gtftk select_by_key -t | wc -l`
      [ "$result" -eq 8 ]
    }
    
     #Keep all transcript with merged intron_size >= 2. Monoexonic are kept.
    @test "select_by_intron_size_10" {
     result=`gtftk get_example -d simple_04 |  gtftk select_by_intron_size -s 2 -a  -m | gtftk select_by_key -t | wc -l`
      [ "$result" -eq 16 ]
    }
   
    '''

    CmdObject(name="select_by_intron_size",
              message="Select transcripts by intron size.",
              parser=make_parser(),
              fun=os.path.abspath(__file__),
              group="selection",
              desc=__doc__,
              notes=__notes__,
              updated=__updated__,
              test=test)

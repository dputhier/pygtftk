from builtins import str

from pybedtools import BedTool

from pygtftk.utils import make_tmp_file
from pygtftk.utils import message


def get_midpoints(self):
    """Returns a bedtools object containing the midpoints of features.


        :Example:
        >>> from pygtftk import bedtool_extension
        >>> fromscratch1 = bedtool_extension.BedTool('chrX 0 100', from_string=True)
        >>> for i in fromscratch1.get_midpoints(): pass
        >>> assert i.start == 49
        >>> assert i.end == 51
        >>> fromscratch1 = bedtool_extension.BedTool('chrX 0 101', from_string=True)
        >>> for i in fromscratch1.get_midpoints(): pass
        >>> assert i.start == 50
        >>> assert i.end == 51
    """

    message("Calling 'get_midpoints'.", type="DEBUG")

    midpoints_bed = make_tmp_file("Midpoints", ".bed")

    n = 1

    for line in self:

        if line.name == ".":
            name = str(n)
        else:
            name = line.name

        if line.strand == ".":
            strand = "+"
        else:
            strand = line.strand

        if line.score == ".":
            score = "."
        else:
            score = line.score

        diff = line.end - line.start

        if diff % 2 != 0:
            # e.g 10-13 (zero based) -> 11-13 one based
            # mipoint is 12 (one-based) -> 11-12 (zero based)
            # e.g 949-1100 (zero based) -> 950-1100 one based
            # mipoint is 1025 (one-based) -> 1024-1025 (zero based)
            # floored division (python 2)...
            line.end = line.start + int(diff // 2) + 1
            line.start = line.end - 1
        else:
            # e.g 10-14 (zero based) -> 11-14 one based
            # mipoint is 12-13 (one-based) -> 11-13 (zero based)
            # e.g 9-5100 (zero based) -> 10-5100 one based
            # mipoint is 2555-2555 (one-based) -> 2554-2555 (zero based)
            # floored division (python 2)...
            # No real center. Take both

            line.start = line.start + int(diff // 2) - 1
            line.end = line.start + 2

        midpoints_bed.write("\t".join([line.chrom,
                                       str(line.start),
                                       str(line.end),
                                       name,
                                       score,
                                       strand]) + "\n")
        n += 1

    midpoints_bed.close()

    return BedTool(fn=midpoints_bed.name)


BedTool.get_midpoints = get_midpoints

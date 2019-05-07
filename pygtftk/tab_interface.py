"""
Class declaration of the TAB object (may be returned by a GTF instance).
"""

import textwrap

import pandas as pd
from cffi import FFI

import pygtftk.utils
from pygtftk.Line import FieldSet
from pygtftk.utils import GTFtkError

ffi = FFI()

# ---------------------------------------------------------------
# The TAB class
# ---------------------------------------------------------------

MAX_REPR_LINE = 6


class TAB(object):
    """A class representation of a tabulated matrix. An end product
    from a GTF object.
    """

    def __init__(self, fn, ptr, dll=None):
        """
        :param ptr: a pointer to a TAB_RESULT structure returned by a GTF instance.
        :param fn: The name of the source GTF file.

        :Example:

        >>> from  pygtftk.utils import get_example_file
        >>> from pygtftk.gtf_interface import GTF
        >>> a_file = get_example_file()
        >>> a_gtf = GTF(a_file)
        >>> a_tab = a_gtf.extract_data("gene_id,start")
        >>> assert a_tab.ncols == 2
        >>> assert len(a_tab) == 70
        >>> assert a_tab.colnames == ['gene_id', 'start']
        >>> assert a_tab.ncols == 2
        >>> assert a_tab.nrows == 70
        >>> assert len(a_tab.as_simple_list(0)) == 70
        >>> assert len(a_tab.as_simple_list(1)) == 70
        >>> assert len(list(set(a_tab.as_simple_list(0)))) == 10
        """

        self.fn = fn
        self._data = ptr
        self._dll = dll

        self.ncols = ptr.nb_columns
        self.nrows = ptr.nb_rows
        self.colnames = []
        for i in range(self.ncols):
            self.colnames += [ffi.string(ffi.cast('char *',
                                                  ptr.column_name[i])).decode()]

    def __iter__(self):
        """
        The iterating function.

        :Example:

        >>> from  pygtftk.utils import get_example_file
        >>> from pygtftk.gtf_interface import GTF
        >>> a_file = get_example_file()[0]
        >>> a_tab = GTF(a_file).extract_data("chrom,gene_id,transcript_id")
        >>> n=0
        >>> for i in a_tab: n += 1
        >>> assert n == 70
        >>> chr = list()
        >>> for i in a_tab: chr += [i[0]]
        >>> assert chr[0] == 'chr1'
        >>> assert chr[::-1][0] == 'chr1'
        """

        for i in range(self.nrows):
            field = FieldSet(self._data.data[i],
                             size=self.ncols,
                             ft_type=self.colnames)
            yield field

    def iter_as_list(self):
        """
        Iterate over the TAB object and return a list of self.nb_columns elements.

        :Example:

        >>> from  pygtftk.utils import get_example_file
        >>> from pygtftk.gtf_interface import GTF
        >>> a_file = get_example_file()
        >>> a_gtf = GTF(a_file)
        >>> a_tab = a_gtf.extract_data("gene_id,start")
        >>> assert next(a_tab.iter_as_list()) == list(a_tab[0])

        """

        for i in range(self.nrows):
            yield [ffi.string(x).decode() for x in self._data.data[i][0:self.ncols]]

    def as_data_frame(self):
        """Convert the TAB object into a dataframe.

        :Example:

        >>> from  pygtftk.utils import get_example_file
        >>> from pygtftk.gtf_interface import GTF
        >>> a_file = get_example_file()
        >>> a_gtf = GTF(a_file)
        >>> a_tab = a_gtf.extract_data("gene_id,start")
        >>> assert a_tab.as_data_frame()["gene_id"].nunique() == 10

        """
        out_list = []
        for i in range(self.nrows):
            out_list += [[ffi.string(x).decode()
                          for x in self._data.data[i][0:self.ncols]]]
        df = pd.DataFrame(out_list, columns=self.colnames)
        return df

    def __getitem__(self, x=None):
        """
        The indexing function.

        :param x: an int.

        :Example:

        >>> from  pygtftk.utils import get_example_file
        >>> from pygtftk.gtf_interface import GTF
        >>> a_file = get_example_file()[0]
        >>> a_tab = GTF(a_file).extract_data("transcript_id,gene_id")
        >>> assert len(a_tab) ==70
        >>> a_tab = GTF(a_file).extract_data("seqid,gene_id")
        >>> assert list(a_tab[len(a_tab)-1]) == ['chr1', 'G0010']
        """

        if not isinstance(x, int):
            raise GTFtkError("This object contains FieldSet objects. "
                             "Index it with integers")
        else:
            if x < self.nrows:
                field = FieldSet(self._data.data[x],
                                 size=self.ncols,
                                 ft_type=self.colnames)
                return field
            else:
                raise GTFtkError("Index out of range.")

    def __repr__(self):
        """Printable representation of the object."""

        if self._data is not None:
            msg = textwrap.dedent("""
            - A TAB object containing {a} lines:
            - File connection: {b}
            """.format(a=self.nrows,
                       b=self.fn))

            n = 0
            for i in self:
                msg += i.format() + "\n"
                n += 1
                if n == MAX_REPR_LINE:
                    break
        else:
            msg = textwrap.dedent("""
            - An empty TAB object.
            - File connection: {b}
            """.format(b=self.fn))

        return msg

    def __len__(self):
        """Object len (number of rows).

        :Example:

        >>> from  pygtftk.utils import get_example_file
        >>> from pygtftk.gtf_interface import GTF
        >>> a_file = get_example_file()[0]
        >>> a_tab = GTF(a_file).extract_data("transcript_id,gene_id")
        >>> assert len(a_tab) == 70

        """

        if self._data is not None:
            return self.nrows
        else:
            return 0

    def write(self, outfile=None, sep='\t'):
        """Write a tab object to a file.

        :Example:

        >>> from  pygtftk.utils import get_example_file
        >>> from pygtftk.gtf_interface import GTF
        >>> from pygtftk.utils import make_tmp_file
        >>> from pygtftk.utils import simple_line_count
        >>> from pygtftk.utils import simple_nb_column
        >>> a_file = get_example_file()[0]
        >>> a_tab = GTF(a_file).extract_data("transcript_id,gene_id")
        >>> out_file = make_tmp_file()
        >>> a_tab.write(out_file)
        >>> out_file.close()
        >>> assert simple_line_count(out_file) == 70
        >>> assert simple_nb_column(out_file) == 2
        """

        for i in self:
            pygtftk.utils.write_properly(sep.join(list(i)), outfile)

    def as_simple_list(self, which_col=0):
        """Convert the selected column of a TAB object into a list.

        :param which_col: The column number.

        :Example:

        >>> from  pygtftk.utils import get_example_file
        >>> from pygtftk.gtf_interface import GTF
        >>> a_file = get_example_file()[0]
        >>> a_gtf = GTF(a_file)[("feature", "gene")]
        >>> a_tab = a_gtf.extract_data("seqid")
        >>> assert list(set(a_tab.as_simple_list(0))) == ['chr1']
        """

        if which_col > self.ncols:
            raise GTFtkError(
                "which_col is greater than the number of contained fields.")

        a_list = list()

        for i in self:
            a_list += [i[which_col]]

        return a_list

    def iterate_with_header(self):
        """Iterate over FieldSets. The first line contains the field names.

        :Example:

        >>> from  pygtftk.utils import get_example_file
        >>> from pygtftk.gtf_interface import GTF
        >>> a_file = get_example_file()[0]
        >>> a_gtf = GTF(a_file).extract_data("all")
        >>> assert next(a_gtf.iterate_with_header())[0] == 'seqid'
        >>> for i in a_gtf.iterate_with_header(): pass
        >>> assert 'CDS_G0010T001' in list(i)
        """
        yield FieldSet(alist=self.colnames)
        for i in range(self.nrows):
            field = FieldSet(self._data.data[i],
                             size=self.ncols,
                             ft_type=self.colnames)
            yield field


if __name__ == "__main__":
    pass

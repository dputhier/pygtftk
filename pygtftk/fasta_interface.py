"""
Class declaration of the FASTA object (may be returned by GTF object methods).
"""

import textwrap

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from cffi import FFI

from pygtftk.Line import FastaSequence

ffi = FFI()


# ---------------------------------------------------------------
# Function def
# ---------------------------------------------------------------


def native_str(x):
    return bytes(x.encode())


# ---------------------------------------------------------------
# The FASTA class
# ---------------------------------------------------------------

MAX_REPR_LINE = 6


class FASTA(object):
    """A class representation of a fasta file connection. A end product
    from a GTF object.
        """

    def __init__(self, fn, ptr, intron, rev_comp):
        """

        :param fn: The file name
        :param ptr: A pointer to the result of GTF._dll.get_sequences()
        :param intron: An int indicating whether intronic regions should be included.
        :param rev_comp: An int indicating whether sequences should be revert-complemented;

        :Example:

        >>> # The __init__ example is provided as doctest
        >>> # Please use get_sequences method from the GTF object.
        >>> from pygtftk.fasta_interface import  FASTA
        >>> from pygtftk.gtf_interface import GTF
        >>> native_str=lambda x:bytes(x.encode())
        >>> from pygtftk.utils import get_example_file
        >>> a_file = get_example_file("simple", "gtf")[0]
        >>> a_gtf = GTF(a_file)
        >>> genome_fa = get_example_file("simple", "fa")[0]
        >>> seq = a_gtf._dll.get_sequences(a_gtf._data, native_str(genome_fa), 0, 0)
        >>> a_fasta = FASTA(a_gtf.fn, seq, False, False)
        >>> assert len(a_fasta) == 15
        """
        self._data = ptr
        self.fn = fn
        self.intron = intron
        self.rev_comp = rev_comp

    def __repr__(self):
        if self._data is not None:
            msg = textwrap.dedent("""
            - A FASTA object containing {a} sequences:
            - File connection: {b}
            - Introns: {c}
            - Reverse-complemented: {d}\n
            """.format(a=str(self._data.nb),
                       b=self.fn,
                       c=self.intron,
                       d=self.rev_comp))

            n = 0
            for i in self:
                msg += i.format() + "\n"
                n += 1
                if n == MAX_REPR_LINE:
                    break
        else:
            msg = textwrap.dedent("""
            - An empty FASTA object.
            - File connection: {b}
            - Introns: {c}
            - Reverse-complemented: {d}\n
            """.format(b=self.fn,
                       c=self.intron,
                       d=self.rev_comp))

        return msg

    def __iter__(self):
        """Iterate over transcript sequences with or without intron depending on arguments provided to GTF.get_sequences().


        :Example:

        >>> from pygtftk.utils import get_example_file
        >>> from pygtftk.gtf_interface import GTF
        >>> a_file = get_example_file("simple", "gtf")[0]
        >>> a_gtf = GTF(a_file)
        >>> genome_fa = get_example_file("simple", "fa")[0]
        >>> a_fa = a_gtf.get_sequences(genome_fa)
        >>> a_dict = dict()
        >>> for i in a_fa: a_dict[i.header] = i.sequence
        """
        for i in range(self._data.nb):
            fs = FastaSequence(self._data.sequence[i])
            yield fs

    def transcript_as_bioseq_records(self):
        """Iterate over transcript sequences with or without intron depending on arguments provided to GTF.get_sequences(). Returns objects of class Bio.SeqRecord.


        :Example:

        >>> from pygtftk.utils import get_example_file
        >>> from pygtftk.gtf_interface import GTF
        >>> a_file = get_example_file("simple", "gtf")[0]
        >>> a_gtf = GTF(a_file)
        >>> genome_fa = get_example_file("simple", "fa")[0]
        >>> a_fa = a_gtf.get_sequences(genome_fa)
        >>> a_list = []
        >>> for i in a_fa.transcript_as_bioseq_records(): a_list += [i]
        >>> assert len(a_list) == 15
        >>> rec = a_list[0]

        """

        description = ''
        if self.intron:
            description += 'immature_transcript'
            if self.rev_comp:
                description += '|rev_comp'
        else:
            description += 'mature_transcript'
            if self.rev_comp:
                description += '|rev_comp'

        for i in range(self._data.nb):
            record = SeqRecord(Seq(ffi.string(self._data.sequence[i].sequence).decode(),
                                   IUPAC.ambiguous_dna),
                               id=ffi.string(
                                   self._data.sequence[i].transcript_id).decode(),
                               name=ffi.string(self._data.sequence[i].gene_id).decode(),
                               description="")
            record.description = description

            yield record

    """
    This should be rewritten carefully.

    def iter_features(self, feat="exon", no_feat_name=False):
        Iterate over features.

        :param feat: The feature type (exon, CDS).
        :param no_feat_name: Don't write the feature type in the header.

        :Example:

        >>> # Load sequence
        >>> from pygtftk.utils import get_example_file
        >>> from pygtftk.gtf_interface import GTF
        >>> a_file = get_example_file("simple", "gtf")[0]
        >>> a_gtf = GTF(a_file)
        >>> genome_fa = get_example_file("simple", "fa")[0]
        >>> a_fa = a_gtf.get_sequences(genome_fa, intron=True)
        >>> # Get exons
        >>> a_list = []
        >>> for i in a_fa.iter_features("exon"): a_list += [i]
        >>> assert len(a_list) == 25
        >>> assert [rec.sequence for rec in a_list if rec.transcript_id == 'G0008T001'] == ['cat', 'gcgct']
        >>> assert [rec.sequence for rec in a_list if rec.transcript_id == 'G0004T001'] == ['atct', 'g', 'gcg']
        >>> # Get exons not reverse-complemented
        >>> a_fa = a_gtf.get_sequences(genome_fa, intron=True, rev_comp=False)
        >>> a_list = []
        >>> for i in a_fa.iter_features("exon"): a_list += [i]
        >>> assert len(a_list) == 25
        >>> assert [rec.sequence for rec in a_list if rec.transcript_id == 'G0008T001'] == ['agcgc', 'atg']
        >>> assert [rec.sequence for rec in a_list if rec.transcript_id == 'G0004T001'] == ['atct', 'g', 'gcg']

   

        if not self.intron:
            GTFtkError("Can't iter_features() if self.intron is False.")

        for i in range(self._data.nb):
            b = self._data.sequence[i]
            seq = ffi.string(b.sequence).decode()
            rg = list(range(b.features.nb))
            for j in rg:
                c = b.features.feature[j]
                name = ffi.string(c.name).decode()
                tx_id = ffi.string(b.transcript_id).decode()
                gn_id = ffi.string(b.gene_id).decode()

                if name == feat:
                    s = c.tr_start
                    e = c.tr_end

                    if no_feat_name:
                        name = [tx_id,
                                gn_id]
                    else:
                        name = [name,
                                tx_id,
                                gn_id]
                    fs = FastaSequence(alist=[">" + "|".join(name),
                                              seq[s:e + 1],
                                              ffi.string(b.seqid).decode(),
                                              b.strand,
                                              gn_id,
                                              tx_id,
                                              c.start,
                                              c.end,
                                              name])
                    yield fs

    def enumerate(self, type="any", add_feature=False):
        Iterate transcript wise and return tx_id and a list of exons
        for i in range(self._data.nb):
            tx = self._data.sequence[i]
            tx_seq = tx.sequence
            seq_list = []
            feat_list = []
            for j in range(tx.features.nb):
                feature = tx.features.feature[j]
                if feature.name == type or type == "any":
                    start = feature.tr_start
                    end = feature.tr_end
                    seq_list += [tx_seq[start:(end + 1)]]
                    if add_feature:
                        feat_list += [feature.name]
            if add_feature:
                yield tx.header, feat_list, seq_list
            else:
                yield tx.header, seq_list

    def as_dict(self, feat="exon"):
        Iterate feature wise

        :param feat: The feature type (exon, CDS).

        :Example:

        >>> from pygtftk.utils import get_example_file
        >>> from pygtftk.gtf_interface import GTF
        >>> a_file = get_example_file("simple", "gtf")[0]
        >>> a_gtf = GTF(a_file)
        >>> genome_fa = get_example_file("simple", "fa")[0]
        >>> a_fa = a_gtf.get_sequences(genome_fa, intron=True)
        >>> a_dict = a_fa.as_dict(feat="exon")
        >>> assert a_dict[('G0001', 'G0001T002', 'chr1', 125, 138, '+', 'exon')] == 'cccccgttacgtag'
        >>> assert a_dict[('G0003', 'G0003T001', 'chr1', 50, 54, '-', 'exon')] == 'gcttg'
        >>> assert a_dict[('G0003', 'G0003T001', 'chr1', 57, 61, '-', 'exon')] == 'aatta'
        >>> assert a_dict[('G0004', 'G0004T001', 'chr1', 71, 71, '+', 'exon')] == 'g'
        >>> assert a_dict[('G0004', 'G0004T001', 'chr1', 65, 68, '+', 'exon')] == 'atct'
        >>> a_dict = a_fa.as_dict(feat="CDS")
        >>> assert a_dict[('G0001', 'G0001T002', 'chr1', 125, 130, '+', 'CDS')] == 'cccccg'
        >>> assert a_dict[('G0003', 'G0003T001', 'chr1', 50, 52, '-', 'CDS')] == 'ttg'
        >>> assert a_dict[('G0004', 'G0004T001', 'chr1', 65, 67, '+', 'CDS')] == 'atc'
        >>> assert a_dict[('G0004', 'G0004T002', 'chr1', 71, 71, '+', 'CDS')] == 'g'
        >>> assert a_dict[('G0004', 'G0004T002', 'chr1', 74, 75, '+', 'CDS')] == 'gc'
        >>> assert a_dict[('G0006', 'G0006T001', 'chr1', 22, 25, '-', 'CDS')] == 'acat'
        >>> assert a_dict[('G0006', 'G0006T001', 'chr1', 28, 30, '-', 'CDS')]  == 'att'
        

        d_out = OrderedDict()

        if not self.intron:
            GTFtkError("Can't iter_features() if self.intron is False.")

        for i in range(self._data.nb):
            b = self._data.sequence[i]
            seq = ffi.string(b.sequence).decode()

            if self.rev_comp:
                if b.strand.decode() == "-":
                    rg = reversed(list(range(b.features.nb)))
                else:
                    rg = list(range(b.features.nb))
            else:
                rg = list(range(b.features.nb))
            for j in rg:
                c = b.features.feature[j]
                name = ffi.string(c.name).decode()
                tx_id = ffi.string(b.transcript_id).decode()
                gn_id = ffi.string(b.gene_id).decode()
                if name == feat:
                    s = c.tr_start
                    e = c.tr_end
                    if b.strand.decode() == "-" and self.rev_comp:
                        end = b.end - s
                        start = b.end - e
                    elif b.strand.decode() == "+":
                        start = b.start + s
                        end = b.start + e
                    elif b.strand.decode() == "":
                        start = b.start + s
                        end = b.start + e
                    else:
                        raise GTFtkError("Strand is undefined")

                    d_out[(gn_id,
                           tx_id,
                           ffi.string(b.seqid).decode(),
                           start,
                           end,
                           b.strand.decode(),
                           name)] = seq[s:e + 1]
        return d_out
    """

    def __getitem__(self, x=None):
        if 0 <= x < self._data.nb:
            return FastaSequence(self._data.data[x])
        else:
            pass

    def __len__(self):

        if self._data is not None:
            return self._data.nb
        else:
            return 0

    def write(self, outputfile):
        """Write the sequences to a fasta file.

        :param outputfile: Output file (file object).

        :Example:

        >>> from pygtftk.utils import get_example_file
        >>> from pygtftk.gtf_interface import GTF
        >>> from pygtftk.utils import make_tmp_file
        >>> from pygtftk.utils import simple_line_count
        >>> a_file = get_example_file("simple", "gtf")[0]
        >>> a_gtf = GTF(a_file)
        >>> genome_fa = get_example_file("simple", "fa")[0]
        >>> a_fa = a_gtf.get_sequences(genome_fa, intron=True)
        >>> tmp_file = make_tmp_file()
        >>> a_fa.write(tmp_file)
        >>> tmp_file.close()
        >>> assert simple_line_count(open(tmp_file.name)) == 30
        """

        for i in range(self._data.nb):
            tx = self._data.sequence[i]
            outputfile.write(ffi.string(tx.header).decode() + "\n")
            outputfile.write(ffi.string(tx.sequence).decode() + "\n")


if __name__ == "__main__":
    pass


Warning about supported GTF file formats
===============================================

.. warning:: Most of the commands of the gtftk suite are designed to handle files in **Ensembl** GTF format and thus require **transcript and gene features/lines** in the GTF. All lines must contain a transcript_id and gene_id value except the **gene feature** that should contain only the gene_id (**see get_example command for an example**). Transcript and gene lines will be used when required to get access to transcript and gene coordinates. This solution was choosen to define a reference GTF file format for (py)gtftk (since Ensembl format is probably the most widely used).

You can use the **convert_ensembl** subcommand to convert your non- (or old) ensembl format to current ensembl format.


Below an example in which we first select only exon features then use **convert_ensembl** to re-generate gene and transcript features using **convert_ensembl** .

.. command-output:: gtftk get_example | gtftk select_by_key -k feature  -v exon | head -n 10
	:shell:


.. command-output:: gtftk get_example | gtftk select_by_key -k feature  -v exon | gtftk  convert_ensembl | head -n 10
	:shell:

**Arguments:**

.. command-output:: gtftk convert_ensembl -h
	:shell:


.. note:: any comment line (*i.e.* starting with #) or empty line in the gtf file will be ignore (discarded) by gtftk.



Naming conventions
----------------------

.. note:: We will use the terms **attribute or key** for any descriptor found in the 9th column (*e.g.* transcript_id) and the term **value** for its associated string (e.g. "NM_334567"). The eight first columns of the GTF file (chrom/seqid, source, type, start, end, score, strand, frame) will be refered as **basic attributes**. In the example below, gene_id is the attribute and 'G0001' is the associated value.

.. command-output:: gtftk get_example| gtftk select_by_key -k feature -v gene| head -1
	:shell:

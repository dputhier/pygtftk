About pygtftk and supported GTF files formats
===============================================

The **Python GTF toolkit (pygtftk) package** is intented to ease handling of GTF (Gene Transfer Format) files. The pygtftk package is compatible with Python 2.7 and relies on **libgtftk**, a library of functions **written in C**.

The package comes with a set of **UNIX commands** that can be accessed through the **gtftk  program**. The gtftk program proposes several atomic tools than can be piped to filter, convert, or extract data from GTF files (...). The gtftk set of Unix commands, can be easily extended using a basic plugin architecture. All these aspects are covered in the help section.

While the gtftk Unix program comes with hundreds of unitary and functional tests, it is still upon  active development and may thus suffer from bugs that remain to be discovered. Feel free to post any problem or required enhancement in the issue section of the github repository.



Supported GTF file formats (you must read this section !)
-----------------------------------------------------------------

.. warning:: The gtftk program is designed to handle files in **ensembl** GTF format. This means that the GTF file provided to gtftk **must contain (for most of the commands) transcript and gene features/lines**. All lines must contain a transcript_id and gene_id value except the **gene feature** that should contain only the gene_id (**see get_example command for an example**). Transcript and gene lines will be used when required to get access to transcript and gene coordinates. This solution was choosen to define a reference GTF file format for gtftk (since Ensembl format is probably the most widely used).

You may use the *convert_ensembl* command to convert your non- (or old) ensembl format to current ensembl format.


Below an example in which we first select only exon features then use *convert_ensembl* to re-generate gene and transcript features.

.. command-output:: gtftk get_example | gtftk select_by_key -k feature  -v exon | head -n 10
	:shell:


.. command-output:: gtftk get_example | gtftk select_by_key -k feature  -v exon | gtftk  convert_ensembl | head -n 10
	:shell:

**Arguments:**

.. command-output:: gtftk convert_ensembl -h
	:shell:


.. note:: any comment line (starting with #) or empty line in the gtf file will be ignore (discarded) by gtftk.



Naming conventions
----------------------

.. note:: We will use the terms **attribute or key** for any descriptor found in the 9th column (*e.g.* transcript_id) and the term **value** for its associated string (e.g. "NM_334567"). The eight first columns of the GTF file (chrom/seqid, source, type, start, end, score, strand, frame) will be refered as **basic attributes**. In the example below, gene_id is the attribute and 'G0001' is the associated value.

.. command-output:: gtftk get_example| gtftk select_by_key -k feature -v gene| head -1
	:shell:


------------------------------------------------------------------------------------------------------------------

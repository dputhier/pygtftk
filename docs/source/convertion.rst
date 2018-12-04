Commands from section 'convertion'
-----------------------------------

convert
~~~~~~~~~~~~~~~~~~~~~~

**Description:** This command can be used to convert to various formats. Currently only a limited number is supported.

* **bed**:  classical bed6 format.
* **bed6**: classical bed6 format.
* **bed3**: bed3 format.


**Example:** Get the gene features and convert them to bed6.

.. command-output:: gtftk get_example | gtftk select_by_key -k feature -v gene | gtftk convert -n gene_id  -f bed6| head -n 3
	:shell:


**Arguments:**

.. command-output:: gtftk convert -h
	:shell:

------------------------------------------------------------------------------------------------------------------

tabulate
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Extract key/values from the GTF and convert them to tabulated format. When requesting coordinates they will be provided in 1-based format.


**Example:** Simply get the list of transcripts and gene.

.. command-output:: gtftk get_example -f gtf | gtftk select_by_key -k feature -v transcript| gtftk tabulate -k gene_id,transcript_id -s "|"
	:shell:

.. warning:: By default tabulate will discard any line for which one of the selected key is not defined. Use -x (--accept-undef) to print them.


**Arguments:**

.. command-output:: gtftk tabulate -h
	:shell:

------------------------------------------------------------------------------------------------------------------


bed_to_gtf
~~~~~~~~~~~~~~~~~~~~~~


**Description:** Convert a bed file to gtf-like format.

**Example:**

.. command-output:: gtftk get_example |gtftk convert| gtftk bed_to_gtf -t transcript | head -n 5
	:shell:


**Arguments:**

.. command-output:: gtftk bed_to_gtf -h
	:shell:


------------------------------------------------------------------------------------------------------------------


convert_ensembl
~~~~~~~~~~~~~~~~~~~~~~


**Description:** Convert the GTF file to ensembl format. Essentially add 'transcript'/'gene' features.

**Example:** Delete gene and transcript feature. Regenerate them.

.. command-output:: gtftk get_example | gtftk select_by_key -k feature -v gene,transcript -n| gtftk convert_ensembl | gtftk select_by_key -k gene_id -v G0001
	:shell:


**Arguments:**

.. command-output:: gtftk bed_to_gtf -h
	:shell:

Commands from section 'selection'
---------------------------------


select_by_key
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Extract lines from the gtf based on key and values.


**Example:** Select some features (genes) then some gene_id.

.. command-output:: gtftk get_example |gtftk select_by_key -k feature -v gene | gtftk select_by_key -k gene_id -v G0002,G0003,G0004
	:shell:


**Arguments:**

.. command-output:: gtftk select_by_key -h
	:shell:

------------------------------------------------------------------------------------------------------------------

select_by_regexp
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Select lines by testing values of a particular key with a regular expression

**Example:** Select lines corresponding to gene_names matching the regular expression 'G.*9$'.

.. command-output:: gtftk get_example |  gtftk select_by_regexp -k gene_id -r 'G.*9$'
	:shell:

**Arguments:**

.. command-output:: gtftk select_by_regexp -h
	:shell:

------------------------------------------------------------------------------------------------------------------

select_by_intron_size
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Delete genes containing an intron whose size is below s. If -m is selected, any gene whose sum of intronic region length is above s is deleted. Monoexonic genes are kept.

**Example:** Some genes having transcripts containing an intron whose size is below 80 nucleotides

.. command-output:: gtftk get_example -d mini_real |  gtftk select_by_intron_size -s 80 -vd | gtftk intron_sizes | gtftk tabulate -k gene_name,transcript_id,intron_sizes -Hun
	:shell:

**Arguments:**

.. command-output:: gtftk select_by_regexp -h
	:shell:

------------------------------------------------------------------------------------------------------------------

select_by_max_exon_nb
~~~~~~~~~~~~~~~~~~~~~~

**Description:** For each gene select the transcript with the highest number of exons.


**Example:** Select lines corresponding to gene_names matching the regular expression 'BCL.*'.

.. command-output:: gtftk get_example |  gtftk select_by_max_exon_nb | gtftk select_by_key -t
	:shell:

**Arguments:**

.. command-output:: gtftk select_by_max_exon_nb -h
	:shell:


------------------------------------------------------------------------------------------------------------------

select_by_loc
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Select transcripts/gene overlapping a given locations. A transcript is defined here as the genomic region from TSS to TTS including introns. This function will return the transcript and all its associated elements (exons, utr,...) even if only a fraction (e.g intron) of the transcript is overlapping the feature. If -/-ft-type is set to 'gene' returns the gene and all its associated elements.

**Example:** Select transcripts at a given location.

.. command-output:: gtftk get_example | gtftk select_by_key -k feature -v transcript | gtftk  select_by_loc -l chr1:10-15
	:shell:

**Arguments:**

.. command-output:: gtftk select_by_loc -h
	:shell:

------------------------------------------------------------------------------------------------------------------

select_by_nb_exon
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Select transcripts based on the number of exons.

**Example:**

.. command-output::  gtftk get_example |  gtftk select_by_nb_exon -m 2 | gtftk nb_exons| gtftk select_by_key -t
	:shell:

**Arguments:**

.. command-output:: gtftk select_by_nb_exon -h
	:shell:


------------------------------------------------------------------------------------------------------------------


select_by_numeric_value
~~~~~~~~~~~~~~~~~~~~~~~~~

**Description:** Select lines from a GTF file based on a boolean test on numeric values.

**Example:**

.. command-output:: gtftk join_attr -i simple.gtf  -j simple.join_mat -k gene_id -m|  gtftk select_by_numeric_value -t 'start < 10 and end > 10 and S1 == 0.5555 and S2 == 0.7' -n ".,?"
	:shell:

**Arguments:**

.. command-output:: gtftk select_by_numeric_value -h
	:shell:


------------------------------------------------------------------------------------------------------------------

random_list
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Select a random list of genes or transcripts.

**Example:** Select randomly 3 transcripts.

.. command-output:: gtftk get_example | gtftk random_list -n 3| gtftk count
	:shell:


**Arguments:**

.. command-output:: gtftk random_list -h
	:shell:

------------------------------------------------------------------------------------------------------------------

random_tx
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Select randomly up to m transcript for each gene.

**Example:** Select randomly 1 transcript per gene (*-m 1*).

.. command-output:: gtftk get_example |  gtftk random_tx -m 1| gtftk select_by_key -k feature -v gene,transcript| gtftk tabulate -k gene_id,transcript_id
	:shell:

**Arguments:**

.. command-output:: gtftk random_tx -h
	:shell:

------------------------------------------------------------------------------------------------------------------

rm_dup_tss
~~~~~~~~~~~~~~~~~~~~~~

**Description:** If several transcripts of a gene share the same tss, select only one.

**Example:** Use rm_dup_tss to select transcripts that will be used for mk_matrix -k 5 (see later).

.. command-output:: gtftk get_example |  gtftk rm_dup_tss| gtftk select_by_key -k feature -v transcript
	:shell:


**Arguments:**

.. command-output:: gtftk rm_dup_tss -h
	:shell:


------------------------------------------------------------------------------------------------------------------

select_by_go
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Select genes from a GTF file using a Gene Ontology ID (e.g GO:0050789).

**Example:** Select genes with transcription factor activity from the GTF. They could be used subsequently to test their epigenetic features (see later).

.. command-output:: # gtftk get_example -d mini_real -f gtf| gtftk select_by_go -s hsapiens | gtftk select_by_key -k feature -v gene | gtftk tabulate -k gene_id,gene_name -Hun | head -6
	:shell:

**Arguments:**

.. command-output:: gtftk select_by_go -h
	:shell:


------------------------------------------------------------------------------------------------------------------

select_by_tx_size
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Select transcript based on their size (i.e size of mature/spliced transcript).

**Example:**

.. command-output:: gtftk get_example | gtftk feature_size -t mature_rna |  gtftk select_by_tx_size -m 14 | gtftk tabulate -n -k gene_id,transcript_id,feat_size
	:shell:


**Arguments:**

.. command-output:: gtftk select_by_tx_size -h
	:shell:

------------------------------------------------------------------------------------------------------------------

select_most_5p_tx
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Select the most 5' transcript of each gene.

**Example:**

.. command-output:: gtftk get_example | gtftk select_most_5p_tx | gtftk select_by_key -k feature -v transcript| gtftk tabulate -k gene_id,transcript_id
	:shell:

**Arguments:**

.. command-output:: gtftk select_most_5p_tx -h
	:shell:

------------------------------------------------------------------------------------------------------------------

short_long
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Get the shortest or longest transcript of each gene

**Example:**

.. command-output:: gtftk get_example | gtftk short_long | gtftk select_by_key -k feature -v transcript| gtftk tabulate -k gene_id,transcript_id
	:shell:

**Arguments:**

.. command-output:: gtftk short_long -h
	:shell:

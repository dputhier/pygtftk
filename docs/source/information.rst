Commands from section 'information'
--------------------------------------


apropos
~~~~~~~~~

**Description:** Search in all command description files those related to a user-defined keyword.

**Example:** Search all commands related to promoters.

.. command-output:: gtftk apropos -k promoter
    :shell:


**Arguments:**

.. command-output:: gtftk apropos -h
    :shell:


------------------------------------------------------------------------------------------------------------------

retrieve
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Retrieves a GTF file from ensembl.

**Example:** List the available GTF files in ensembl FTP. Bacteria are not listed at the moment.

.. command-output:: # gtftk retrieve -l | head -5
    :shell:

**Example:** Perform basic statistics on Vicugna pacos genomic annotations.

.. command-output:: # gtftk retrieve -s vicugna_pacos -c  -d | gtftk  count -t vicugna_pacos
    :shell:

**Arguments:**

.. command-output:: gtftk retrieve -h
    :shell:


------------------------------------------------------------------------------------------------------------------

get_example
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Get an example GTF file (or any other kind of example available in the installation directory). This command is only provided for demonstration purpose.

We can see from the example below that this gtf file **follows the ensembl format** and contains the **transcript and gene features** (column 3).


**Example:** The very basic (and artificial) example.

.. command-output:: gtftk get_example| head -2
    :shell:


let's get all files from the *simple* dataset.

.. command-output:: gtftk get_example -d simple -f '*'
    :shell:

**Arguments:**

.. command-output:: gtftk get_example -h
    :shell:

------------------------------------------------------------------------------------------------------------------

add_exon_nb
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Add exon number transcript-wise (based on 5' to 3' orientation).

**Example:**

.. command-output:: gtftk get_example | gtftk add_exon_nb  -k exon_number | gtftk select_by_key -k feature -v exon | gtftk tabulate -k chrom,start,end,exon_number,transcript_id | head -n 20
    :shell:


**Arguments:**


.. command-output:: gtftk add_exon_nb -h
    :shell:


------------------------------------------------------------------------------------------------------------------

count
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Count the number of features (transcripts, genes, exons, introns).

**Example:**

.. command-output:: gtftk  get_example | gtftk count  -t example_gtf
    :shell:


**Arguments:**

.. command-output:: gtftk count -h


------------------------------------------------------------------------------------------------------------------

count_key_values
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Count the number of values for a set of keys.


**Example:** Count the number of non-redondant entries for chromosomes and transcript_id.


.. command-output:: gtftk get_example | gtftk count_key_values -k chrom,transcript_id -u
    :shell:


**Arguments:**

.. command-output:: gtftk count_key_values -h


------------------------------------------------------------------------------------------------------------------

get_attr_list
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Get the list of attributes from a GTF file.

**Example:** Get the list of attributes in the "simple" dataset.

.. command-output:: gtftk get_example | gtftk get_attr_list
    :shell:

**Arguments:**

.. command-output:: gtftk get_attr_list -h


------------------------------------------------------------------------------------------------------------------

get_attr_value_list
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Get the list of values observed for an attributes.


**Example:** Get the number of time each gene_id is used.

.. command-output:: gtftk get_example | gtftk get_attr_value_list -k gene_id -c -s ';'
    :shell:


**Arguments:**

.. command-output:: gtftk get_attr_value_list -h


------------------------------------------------------------------------------------------------------------------

get_feature_list
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Get the list of features enclosed in the GTF.

**Example:** Get the list of features enclosed in the GTF.

.. command-output:: gtftk get_example | gtftk get_feature_list
    :shell:

**Arguments:**

.. command-output:: gtftk get_feature_list -h


------------------------------------------------------------------------------------------------------------------

nb_exons
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Count the number of exons and add it as a novel key/value. Output may also be in text format if requested.

**Example:**

.. command-output:: gtftk  get_example -f gtf | gtftk nb_exons | head -n 5
    :shell:

**Arguments:**

.. command-output:: gtftk nb_exons -h
    :shell:


------------------------------------------------------------------------------------------------------------------

nb_transcripts
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Count the number of transcript per gene.

**Example:** Count the number of transcript per gene.

.. command-output:: gtftk get_example |  gtftk nb_transcripts  | gtftk select_by_key -g
    :shell:


**Arguments:**

.. command-output:: gtftk nb_transcripts -h
    :shell:

------------------------------------------------------------------------------------------------------------------

seqid_list
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Return the chromosome list.

**Example:** Return the chromosome list.

.. command-output:: gtftk get_example |  gtftk seqid_list
    :shell:


**Arguments:**

.. command-output:: gtftk seqid_list -h
    :shell:

------------------------------------------------------------------------------------------------------------------

tss_dist
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Compute the distance between TSSs of pairs of gene transcripts. The tss_num_1 and tss_num_2 columns contains the TSSs number (for transcript_id_1 and transcript_id_2 respectively). Numering starts from 1 (most 5' TSS  for a gene) to the number of different TSS coordinates. Two or more transcripts will have the same tss_num if they share a TSS.

**Example:** An example on the mini_real dataset.

.. command-output:: gtftk get_example -d mini_real |  gtftk tss_dist | head -n 10
    :shell:


**Arguments:**

.. command-output:: gtftk tss_dist -h
    :shell:

------------------------------------------------------------------------------------------------------------------


feature_size
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Get the size and limits (start/end) of features enclosed in the GTF. If bed format is requested returns the limits in bed format and the size as a score. Otherwise output GTF file with 'feat_size' as a new key and size as value


**Example:** Add trancript size (mature RNA) to the gtf.

.. command-output:: gtftk get_example | gtftk feature_size -t mature_rna | gtftk select_by_key -k feature -v transcript | head -n 5
    :shell:

**Arguments:**

.. command-output:: gtftk feature_size -h
    :shell:



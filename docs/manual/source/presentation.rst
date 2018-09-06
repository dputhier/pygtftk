Help on gtftk Unix commands
============================



Main parser arguments of gtftk
-------------------------------


Getting help with -h
~~~~~~~~~~~~~~~~~~~~~

The -h argument can be used to get a synopsis for implemented commands.

.. command-output:: gtftk -h
	:shell:


------------------------------------------------------------------------------------------------------------------


Activating Bash completion
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The code provided below can be useful to activate bash completion.

.. code-block:: bash

  # Use the -b argument of gtftk
  # This will produce a script that you
  # should store in your .bashrc
  gtftk -b

Or alternatively

.. code-block:: bash

  echo "" >> ~/.bashrc
  gtftk -b >> ~/.bashrc 



------------------------------------------------------------------------------------------------------------------


Getting the list of funtional tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One can access to the list of implemented  tests through the -p/--plugin-tests arguments. These tests may be run using `bats <https://github.com/sstephenson/bats>`_ (Bash Automated Testing System).


.. code-block:: bash

  # gtftk --plugin-tests



------------------------------------------------------------------------------------------------------------------


Command-wide arguments
--------------------------

**Description:** The following arguments are available in almost all gtftk commands :

- -h, --help : Refers to argument list and details.
- -i, --inputfile: Refers to the input file (may be <stdin>).
- -o, --outputfile: Refers to the output file (may be <stdout>).
- -D, --no-date: Do not add date to output file names.
- -C, --add-chr: Add 'chr' to chromosome names before printing output.
- -V, --verbosity: Increase output verbosity (can take value from 0 to 4).
- -K --tmp-dir: Keep all temporary files into this folder. 
- -L, --logger-file: Store the values of all command line arguments into a file.


------------------------------------------------------------------------------------------------------------------

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

**Description:** Retrieve a GTF file from ensembl.

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


**Example:** The very basic (and artificial example).

.. command-output:: gtftk get_example| head -2
	:shell:


**Example:** A more realistic example containing a subset of transcript (n=8531) corresponding to 1058 genes from human annotation. 

.. command-output:: gtftk get_example -d mini_real | gtftk count
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

.. command-output:: gtftk  get_example -f gtf | gtftk add_exon_nb  | gtftk select_by_key -k feature -v exon
	:shell:

.. command-output:: gtftk get_example -f gtf | gtftk add_exon_nb  -k exon_number | gtftk select_by_key -k feature -v exon | gtftk tabulate -k chrom,start,end,exon_number,transcript_id | head -n 20
	:shell:

**Arguments:**

.. command-output:: gtftk add_exon_nb -h 
	:shell:


------------------------------------------------------------------------------------------------------------------

count
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Count the number of features (transcripts, genes, exons, introns).

**Example:**

.. command-output:: gtftk  get_example -f gtf | gtftk count  -t example_gtf
	:shell:


**Arguments:**

.. command-output:: gtftk count -h


------------------------------------------------------------------------------------------------------------------

count_key_values
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Count the number values for a set of keys.

**Example:** Count the number of time gene_id and transcript_id appear in the GTF file.

.. command-output:: gtftk get_example | gtftk count_key_values -k gene_id,transcript_id
	:shell:

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

**Example:** Get the list of values observed for transcript_id.

.. command-output:: gtftk get_example | gtftk get_attr_value_list -k transcript_id
	:shell:

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

.. command-output:: gtftk  get_example -f gtf | gtftk nb_exons  | gtftk select_by_key -k feature -v transcript | head -n 5
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

**Description:** Returns the chromosome list.

**Example:** Returns the chromosome list.

.. command-output:: gtftk get_example |  gtftk seqid_list
	:shell:


**Arguments:**

.. command-output:: gtftk seqid_list -h
	:shell:

------------------------------------------------------------------------------------------------------------------

tss_dist
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Computes the distance between TSSs of pairs of gene transcripts. The tss_num_1/tss_num_1 columns contains the numbering of TSSs (transcript_id_1 and transcript_id_2 respectively) for each gene. Numering starts from 1 (most 5' TSS) to the number of different TSS coordinates. Two or more transcripts will have the same tss_num if they share a TSS.

**Example:** Returns the chromosome list.

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

**Example:** Add trancript size (genomic coverage) to the gtf.

.. command-output:: gtftk get_example | gtftk feature_size -t transcript | gtftk select_by_key -k feature -v transcript | head -n 5
	:shell:

**Example:** Get exon size and limits in BED format.

.. command-output:: gtftk get_example | gtftk feature_size  -t exon -b -n feature,exon_id,gene_id| head -n 5
	:shell:


**Arguments:**

.. command-output:: gtftk feature_size -h
	:shell:


------------------------------------------------------------------------------------------------------------------

Commands from section 'Editing'
----------------------------------


add_prefix
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Add a prefix (or suffix) to one of the attribute value (*e.g.* gene_id)

**Example:**

.. command-output:: gtftk get_example| gtftk add_prefix -k transcript_id -t "novel_"| head -2
	:shell:

.. command-output:: gtftk get_example| gtftk add_prefix -k transcript_id -t "_novel" -s | head -2
	:shell:

**Arguments:**

.. command-output:: gtftk add_prefix -h
	:shell:

------------------------------------------------------------------------------------------------------------------

del_attr
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Delete an attribute and its corresponding values.

**Example:**

.. command-output:: gtftk get_example | gtftk del_attr -k transcript_id,gene_id,exon_id | head -3
	:shell:

.. command-output:: gtftk get_example | gtftk del_attr -v  -k transcript_id,gene_id | head -3 # delete all but transcript_id and gene_id
	:shell:

**Arguments:**

.. command-output:: gtftk del_attr -h
	:shell:


------------------------------------------------------------------------------------------------------------------

join_attr
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Add attributes from a file. This command can be used to import additional key/values into the gtf (e.g CPAT for coding potential, DESeq for differential analysis,...). The imported file can be in 2 formats (2 columns or matrix):

- With a 2-columns file:

  - value for joining (transcript_id or gene_id or ...).
  - corresponding value.

- With a matrix (see -m):

  - rows corresponding to joining keys (transcript_id or gene_id or...).
  - columns corresponding to novel attributes name.
  - Each cell of the matrix is a value for the corresponding attribute.


**Example:** With a 2-columns file.

.. command-output:: gtftk get_example -f join > simple_join.txt
	:shell:

.. command-output:: cat simple_join.txt
	:shell:

.. command-output::  gtftk get_example -f gtf | gtftk join_attr -k gene_id -j simple_join.txt -n a_score -t gene| gtftk select_by_key -k feature -v gene
	:shell:

**Example:** With a matrix

.. command-output:: gtftk get_example -f join_mat  >  simple_join_mat.txt
	:shell:

.. command-output:: cat simple_join_mat.txt
	:shell:

.. command-output:: gtftk get_example -f gtf | gtftk join_attr -k gene_id -j simple_join_mat.txt -m -t gene| gtftk select_by_key -k feature -v gene
	:shell:


**Arguments:**

.. command-output:: gtftk join_attr -h
	:shell:


------------------------------------------------------------------------------------------------------------------

join_multi_file
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Join attributes from mutiple files.


**Example:** Add key/value to gene feature.

.. command-output:: gtftk get_example |  gtftk join_multi_file -k gene_id -t gene simple.join_mat_2 simple.join_mat_3| gtftk select_by_key -g
	:shell:

**Arguments:**

.. command-output:: gtftk join_multi_file -h
	:shell:



------------------------------------------------------------------------------------------------------------------

merge_attr
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Merge a set of attributes into a destination attribute.


**Example:** Merge gene_id and transcript_id into a new key associated to transcript features.

.. command-output:: gtftk get_example |  gtftk merge_attr -k transcript_id,gene_id -d txgn_id -s "|" -f transcript | gtftk select_by_key -t
	:shell:

**Example:** Merge gene_id and transcript_id into an existing key (transcript_id) associated to transcript features.

.. command-output:: gtftk get_example |  gtftk merge_attr -k transcript_id,gene_id -d transcript_id -s "|" -f transcript | gtftk tabulate -k feature,transcript_id -x | head -n 6
	:shell:

**Example:** Merge gene_id and transcript_id into an existing key (transcript_id) associated to all features.

.. command-output:: gtftk get_example |  gtftk merge_attr -k transcript_id,gene_id -d transcript_id -s "|" -f "*" | gtftk tabulate -k feature,transcript_id -x | head -n 6
	:shell:


**Arguments:**

.. command-output:: gtftk join_multi_file -h
	:shell:


------------------------------------------------------------------------------------------------------------------


discretize_key
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Create a new key by discretizing a numeric key. This can be helpful to create new classes on the fly that can be used subsequently.
The default is to create equally spaced interval. The intervals can also be created by computing the percentiles (-p).


**Example:** Let say we have the following matrix giving expression level of genes (rows) in samples (columns). We could join this information to the GTF and later choose to transform key *S1* into a new discretized key *S1_d*. We may apply particular labels to this factor using *-l*.


.. command-output:: gtftk get_example |  gtftk join_attr -j simple.join_mat -k gene_id -m | gtftk discretize_key -k S1 -d S1_d -n 2 | gtftk select_by_key -k feature -v gene
	:shell:



.. command-output:: gtftk get_example |  gtftk join_attr -j simple.join_mat -k gene_id -m | gtftk discretize_key -k S1 -d S1_d -n 2 -l A,B  | gtftk select_by_key -k feature -v gene
	:shell:

**Arguments:**

.. command-output:: gtftk discretize_key -h
	:shell:

------------------------------------------------------------------------------------------------------------------

Commands from section 'selection'
---------------------------------


select_by_key
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Extract lines from the gtf based on key and values.


**Example:** Select some features (genes) then some gene_id.

.. command-output:: gtftk get_example |gtftk select_by_key -k feature -v gene | gtftk select_by_key -k gene_id -v G0002,G0003,G0004
	:shell:

**Example:** Select gene list in column 1 of file simple_join.txt.

.. command-output:: gtftk get_example -f join > simple_join.txt ; gtftk get_example| gtftk select_by_key -f simple_join.txt -c 1 -k gene_id | gtftk tabulate -k gene_id -Hun
	:shell:

**Example:** Select the gene list enclosed in column 1 of file simple_join.txt. Ask for bed format.

.. command-output:: gtftk get_example -f join > simple_join.txt ; gtftk get_example| gtftk select_by_key -f simple_join.txt -c 1 -k gene_id -b
	:shell:

**Example:** Select all but genes in column 1 of file simple_join.txt.

.. command-output:: gtftk get_example -f join > simple_join.txt ; gtftk get_example| gtftk select_by_key -f simple_join.txt -c 1 -k gene_id -n | gtftk tabulate -k gene_id -Hun
	:shell:

**Arguments:**

.. command-output:: gtftk select_by_key -h
	:shell:

------------------------------------------------------------------------------------------------------------------

select_by_regexp
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Select lines by testing values of a particular key with a regular expression

**Example:** Select lines corresponding to gene_names matching the regular expression 'BCL.*'.

.. command-output:: gtftk get_example -d mini_real |  gtftk select_by_regexp -k gene_name -r "BCL.*" | gtftk tabulate -Hun -k gene_name
	:shell:

**Arguments:**

.. command-output:: gtftk select_by_regexp -h
	:shell:

------------------------------------------------------------------------------------------------------------------

select_by_intron_size
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Delete genes containing an intron whose size is below s. If -m is selected, any gene whose sum of intronic region length is above s is deleted. Monoexonic genes are kept.

**Example:** Select lines corresponding to gene_names matching the regular expression 'BCL.*'.

.. command-output:: gtftk get_example -d mini_real |  gtftk select_by_regexp -k gene_name -r "BCL.*"  | gtftk tabulate -Hun -k gene_name
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

.. command-output:: gtftk get_example -d mini_real -f gtf| gtftk select_by_go -s hsapiens | gtftk select_by_key -k feature -v gene | gtftk tabulate -k gene_id,gene_name -Hun | head -6
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

.. command-output:: gtftk get_example | gtftk feature_size -t mature_rna |  gtftk select_by_tx_size -m 11 | gtftk tabulate -n -k gene_id,transcript_id,feat_size
	:shell:

.. command-output:: gtftk get_example -d mini_real | gtftk feature_size -t mature_rna |  gtftk select_by_tx_size -m 8000  -M 1000000000 | gtftk tabulate -n -k gene_id,transcript_id,feat_size -H  | sort -k3,3n | tail -n 10
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

------------------------------------------------------------------------------------------------------------------



Commands from section 'convertion'
-----------------------------------

convert
~~~~~~~~~~~~~~~~~~~~~~

**Description:** This command can be used to convert to various formats. Currently only a limited number is supported.

* **bed**:  classical bed6 format.
* **bed6**: classical bed6 format.
* **bed3**: bed3 format.


**Example:** Get the gene features and convert them to bed6.

.. command-output:: gtftk get_example | gtftk select_by_key -k feature -v gene | gtftk convert -n gene_id | head -n 3
	:shell:

**Example:** Get the gene features and convert them to bed3.

.. command-output:: gtftk get_example | gtftk select_by_key -k feature -v gene | gtftk convert -f bed3 | head -n 3
	:shell:

**Example:** Get the exonic features and convert them to bed3.

.. command-output:: gtftk get_example | gtftk select_by_key -k feature -v exon | gtftk convert -n gene_id,transcript_id,exon_id | head -3
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


**Example:** Join novel attributes (see **join_attr examples**) and convert the resulting GTF stream to tab format

.. command-output:: gtftk get_example -f gtf | gtftk join_attr -k gene_id -j simple_join.txt -n a_score -t gene| gtftk select_by_key -k feature -v gene| gtftk tabulate -k feature,start,end,seqid,gene_id,a_score
	:shell:

**Example:** You may also delete the header, ask for non redondant lines and delete any lines containing not-available values ('.').

.. command-output:: gtftk get_example -f gtf | gtftk join_attr -k gene_id -j simple_join.txt -n a_score -t gene| gtftk select_by_key -k feature -v gene| gtftk tabulate -k feature,start,end,seqid,gene_id,a_score -Hun
	:shell:


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


------------------------------------------------------------------------------------------------------------------


Commands from section 'annotation'
------------------------------------


closest_genes
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Find the n closest genes for each transcript.

**Example:**

.. command-output:: gtftk get_example |  bedtools sort | gtftk closest_genes -f
	:shell:


**Arguments:**

.. command-output:: gtftk closest_genes -h
	:shell:


overlapping
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Find transcripts whose body/TSS/TTS region extended in 5' and 3' (-u/-d) overlaps with any transcript from another gene. Strandness is not considered by default. Used --invert-match to find those that do not overlap. If --annotate-gtf is used, all lines of the input GTF file will be printed and a new key containing the list of overlapping transcripts will be added to the transcript features/lines (key will be 'overlapping_*' with * one of body/TSS/TTS). The --annotate-gtf and --invert-match arguments are mutually exclusive.


**Example:** Find transcript whose promoter overlap transcript from other genes.

.. command-output:: gtftk get_example -f chromInfo > simple_join_chromInfo.txt;  gtftk get_example | gtftk overlapping -c simple_join_chromInfo.txt -t promoter -u 10 -d 10 -a    | gtftk select_by_key -k feature -v transcript | gtftk tabulate -k transcript_id,overlap_promoter_u0.01k_d0.01k | head
	:shell:

**Example:** Find transcript whose tts overlap transcript from other genes (on the other strand).


.. command-output:: gtftk get_example -f chromInfo > simple_join_chromInfo.txt;  gtftk get_example | gtftk overlapping -c simple_join_chromInfo.txt -t tts -u 30 -d 30 -a -S     | gtftk select_by_key -k feature -v transcript | gtftk tabulate -k transcript_id,overlap_tts_u0.03k_d0.03k | head
	:shell:

**Arguments:**

.. command-output:: gtftk overlapping -h
	:shell:

------------------------------------------------------------------------------------------------------------------

divergent
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Find transcript with divergent promoters. These transcripts will be defined here
as those whose promoter region (defined by -u/-d) overlaps with the tss of
another gene in reverse/antisens orientation. This may be useful to select
coding genes in head-to-head orientation or LUAT as described in "Divergent
transcription is associated with promoters of transcriptional regulators"
(Lepoivre C, BMC Genomics, 2013). The ouput is a GTF with an additional key
('divergent') whose value is set to '.' if the gene has no antisens transcript
in its promoter region. If the gene has an antisens transcript in its promoter
region the 'divergent' key is set to the identifier of the transcript whose tss
is the closest relative to the considered promoter. The tss to tss distance is
also provided as an additional key (dist_to_divergent).


**Example:** Flag divergent transcripts in the example dataset. Select them and produce a tabulated output.

.. command-output:: gtftk get_example -f chromInfo > simple_join_chromInfo.txt;  gtftk get_example |  gtftk divergent -c simple_join_chromInfo.txt -u 10 -d 10| gtftk select_by_key -k feature -v transcript | gtftk tabulate -k transcript_id,divergent,dist_to_divergent | head  -n 7
	:shell:

**Arguments:**

.. command-output:: gtftk divergent -h
	:shell:

------------------------------------------------------------------------------------------------------------------

convergent
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Find transcript with convergent tts. These transcripts will be defined here
as those whose tts region (defined by -u/-d) overlaps with the tts of
another gene in reverse/antisens orientation. The ouput is a GTF with an
additional key ('convergent') whose value is set to '.' if the gene has no
convergent transcript in its tts region. If the gene has an antisens transcript
in its tts region the 'convergent' key is set to the identifier of the
transcript whose tts is the closest relative to the considered tts.
The tts to tts distance is also provided as an additional key (dist_to_convergent).


**Example:** Flag divergent transcripts in the example dataset. Select them and produce a tabulated output.

.. command-output:: gtftk get_example -f chromInfo > simple_join_chromInfo.txt;  gtftk get_example |  gtftk convergent -c simple_join_chromInfo.txt -u 25 -d 25| gtftk select_by_key -k feature -v transcript | gtftk tabulate -k transcript_id,convergent,dist_to_convergent| head -n 4
	:shell:

**Arguments:**

.. command-output:: gtftk convergent -h
	:shell:

------------------------------------------------------------------------------------------------------------------

exon_sizes
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Add a new key to transcript features containing a comma separated list of exon sizes.


**Example:**

.. command-output:: gtftk get_example | gtftk exon_sizes | gtftk select_by_key -t
	:shell:

**Arguments:**

.. command-output:: gtftk exon_sizes -h
	:shell:

------------------------------------------------------------------------------------------------------------------


intron_sizes
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Add a new key to transcript features containing a comma separated list of intron sizes.


**Example:**

.. command-output:: gtftk get_example | gtftk intron_sizes | gtftk select_by_key -t
	:shell:

**Arguments:**

.. command-output:: gtftk intron_sizes -h
	:shell:

------------------------------------------------------------------------------------------------------------------


Commands from section 'coordinates'
-----------------------------------

midpoints
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Get the genomic midpoint of each features: genes, transcripts, exons or introns. Output is currently in bed format only.


**Example:** Get mipoints of all transcripts and exons.

.. command-output:: gtftk get_example | gtftk midpoints -t transcript,exon -n transcript_id,feature | head -n 5
	:shell:


**Arguments:**

.. command-output:: gtftk midpoints -h
	:shell:

------------------------------------------------------------------------------------------------------------------

5p_3p_coord
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Get the 5p or 3p coordinates for each feature (e.g TSS or TTS for a transcript).
Output is bed format.

**Example:** Get the 5p ends of transcripts and exons.

.. command-output:: gtftk get_example | gtftk 5p_3p_coord -t transcript,exon -n transcript_id,gene_id,feature | head -n 5
	:shell:

**Example:** Get the 3p ends of transcripts and exons.

.. command-output:: gtftk get_example | gtftk 5p_3p_coord -t transcript,exon -n transcript_id,gene_id,feature -v -s "^"| head -n 5
	:shell:


**Arguments:**

.. command-output:: gtftk 5p_3p_coord -h
	:shell:

------------------------------------------------------------------------------------------------------------------


intergenic
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Extract intergenic regions. This command requires a chromInfo file to compute
the bed file boundaries. The command will print the coordinates of genomic
regions without transcript features.


**Example:** Simply get intergenic regions.

.. command-output::  gtftk get_example -f chromInfo > simple_join_chromInfo.txt; gtftk get_example |  gtftk intergenic   -c simple_join_chromInfo.txt
	:shell:

**Arguments:**

.. command-output:: gtftk intergenic -h
	:shell:

------------------------------------------------------------------------------------------------------------------

intronic
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Returns a bed file containing the intronic regions. If by_transcript is false
(default), returns merged genic regions with no exonic overlap ("strict" mode).
Otherwise, the intronic regions corresponding to each transcript are returned
(may contain exonic overlap and redundancy).

**Example:** Simply get intronic regions.

.. command-output:: gtftk get_example |  gtftk intronic | head -n 5
	:shell:

**Example:** Intronic regions of each transcript.

.. command-output:: gtftk get_example |  gtftk intronic -b
	:shell:

**Arguments:**

.. command-output:: gtftk intronic -h
	:shell:

------------------------------------------------------------------------------------------------------------------


splicing_site
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Compute the locations of donor and acceptor splice sites. This command will return a single position which corresponds to the most 5' and/or the most 3' intronic region. If the gtf file does not contain exon numbering you can compute it using the
add_exon_nb command. The score column of the bed file contain the number of the closest exon relative to the splice site.

**Example:**

.. command-output:: gtftk get_example | gtftk add_exon_nb -k exon_nbr | gtftk splicing_site  -k exon_nbr| head
	:shell:

**Arguments:**

.. command-output:: gtftk splicing_site -h
	:shell:

------------------------------------------------------------------------------------------------------------------

shift
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Shift coordinates in 3' or 5' direction.

**Example:**

.. command-output:: gtftk get_example|  head -n 1
	:shell:

.. command-output:: gtftk get_example -f chromInfo > simple.chromInfo; gtftk get_example |  gtftk shift -s -10 -c simple.chromInfo | head -n 1
	:shell:


**Arguments:**

.. command-output:: gtftk shift -h
	:shell:


------------------------------------------------------------------------------------------------------------------

Commands from section 'sequence'
---------------------------------


get_tx_seq
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Get transcript sequences in fasta format.

**Example:** Get sequences of transcripts in 5' to 3' orientation

.. command-output:: gtftk get_example -f fa > simple.fa; gtftk get_example | gtftk get_tx_seq -g simple.fa | head -n 4
	:shell:

Note that the format is rather flexible and any combination of key can be exported to the header.

.. command-output:: gtftk get_example | gtftk get_tx_seq -g simple.fa  -l gene_id,transcript_id,feature,chrom,start,end,strand  | head -n 2
	:shell:

You can ask to add explicitly (-e) the name of the keys in the header. Here we also add the size of the mature transcript and the number of exons.

.. command-output:: gtftk get_example | gtftk feature_size -t mature_rna | gtftk nb_exons| gtftk get_tx_seq -g simple.fa -l feature,transcript_id,gene_id,seqid,start,end,feat_size,nb_exons -e | head -n 2
	:shell:

You may use wildcard (path enclosed within quotes) in case the genome is splitted in several chromosome files:

.. command-output:: gtftk get_example |  gtftk get_tx_seq -g '*.fa' -l gene_id,transcript_id,feature,chrom,start,end,strand -s "," | head -n 2
	:shell:

A particular header format that should be compliant with sleuth is also proposed.

.. command-output:: gtftk get_example |  gtftk get_tx_seq -g '*.fa'  -f -n  | head -n 2
	:shell:

**Arguments:**

.. command-output:: gtftk get_tx_seq -h
	:shell:

------------------------------------------------------------------------------------------------------------------

get_feat_seq
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Get feature sequence (e.g exon, UTR...).


**Example:**

.. command-output:: gtftk get_feat_seq -i simple.gtf -g simple.fa  -l feature,transcript_id,start -t  exon -n | head -10
	:shell:

**Arguments:**

.. command-output:: gtftk get_feat_seq -h
	:shell:


------------------------------------------------------------------------------------------------------------------


Commands from section 'coverage'
----------------------------------

coverage
~~~~~~~~

**Description:** Takes a GTF as input to compute bigwig coverage in regions of interest (promoter, transcript body, intron, intron_by_tx, tts...) or a BED6 to focus on user-defined regions. If --n-highest is used the program will compute the coverage of each bigwig based on the average value of the n windows (--nb-window) with the highest coverage values.
Regions were signal can be computed (if GTF file as input) are promoter, tts, introns, intergenic regions or any feature available in the GTF file (transcript, exon, gene...).
If --matrix-out is selected, the signal for each bigwig will be provided in a dedicated column. Otherwise, signal for each bigwig is provided through a dedicated line.


 **Example:**

We will first request a lightweight example dataset.


.. command-output:: gtftk get_example -d mini_real -f '*'
	:shell:


Although we could work on the full dataset, we will focus on transcripts whose promoter region do not overlaps with any transcript from another gene.


.. command-output:: gtftk overlapping -i mini_real.gtf.gz -c hg38.genome  -n > mini_real_noov.gtf
	:shell:


We will select a representative transcript for each gene. Here we will perform this step using random_tx although another interesting choice would be rm_dup_tss.

.. command-output:: gtftk random_tx -i mini_real_noov.gtf  -m 1 -s 123 > mini_real_noov_rnd_tx.gtf
	:shell:

Now we will compute coverage of promoters regions using 3 bigWig files as input.


.. command-output:: gtftk coverage -l H3K4me3,H3K79me2,H3K36me3 -u 5000 -d 5000 -i mini_real_noov_rnd_tx.gtf -c hg38.genome -m transcript_id,gene_name -x ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw > coverage.bed
	:shell:


Now we can have a look at the result:

.. command-output:: head -n 10 coverage.bed
	:shell:


**Arguments:**

.. command-output::  gtftk coverage -h
	:shell:

------------------------------------------------------------------------------------------------------------------


mk_matrix
~~~~~~~~~~

Description: Gtftk implements commands that can be used to produce coverage profiles around genomic features or inside user-defined regions. A coverage matrix need first to  be produced from a bwig using the mk_matrix command.

**Example:**

We will used the same dataset (mini_real.gtf) as produced above (see help on *coverage* command).

We can now create a coverage matrix around TSS/TTS or along the full transcript
(with or without 5' and 3' regions). Provide a BED file as *---inputfile* if you
want to use your own, user-specific, regions.
Will will create tree example datasets:

First we will create a coverage matrix around promoter based on a subset of randomly choose transcripts (one per gene) from the 'mini_real' dataset (see section on the *coverage* command to get info about the construction of the *mini_real_noov_rnd_tx.gtf.gz* dataset).

.. command-output::  gtftk get_example -f '*' -d mini_real
	:shell:


.. command-output:: gtftk mk_matrix -k 5 -i mini_real_noov_rnd_tx.gtf.gz -d 5000 -u 5000 -w 200 -c hg38.genome  -l  H3K4me3,H3K79me,H3K36me3 ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_promoter
	:shell:


Then we will also compute coverage profil around around tts.


.. command-output:: gtftk mk_matrix -k 5 -i mini_real_noov_rnd_tx.gtf.gz -t tts  -d 5000 -u 5000 -w 200 -c hg38.genome  -l  H3K4me3,H3K79me,H3K36me3 ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_tts
	:shell:


The following command compute coverage profil along the whole transcript

.. command-output:: gtftk mk_matrix -k 5 -i mini_real_noov_rnd_tx.gtf.gz -t transcript  -d 5000 -u 5000 -w 200 -c hg38.genome  -l  H3K4me3,H3K79me,H3K36me3 ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_tx
	:shell:


Along the whole transcript but increasing the number of windows dedicated to upstream and downstream regions.


.. command-output:: gtftk mk_matrix -k 5 --bin-around-frac 0.5 -i mini_real_noov_rnd_tx.gtf.gz -t transcript  -d 5000 -u 5000 -w 200 -c hg38.genome  -l  H3K4me3,H3K79me,H3K36me3 ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_tx_2
	:shell:

Along a user defined set of regions (in bed6 format). Here we will used the transcript coordinates in bed format as an example.

.. command-output:: gtftk select_by_key -i mini_real_noov_rnd_tx.gtf.gz -k feature -v transcript | gtftk convert -f bed6 > mini_real_rnd_tx.bed
	:shell:

.. command-output:: gtftk mk_matrix -k 5 --bin-around-frac 0.5 -i mini_real_rnd_tx.bed -t user_regions  -d 5000 -u 5000 -w 200 -c hg38.genome  -l  H3K4me3,H3K79me,H3K36me3 ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_user_def
	:shell:

And finally using a set of single nucleotides coordinates that will be extend (-u/-d) and assessed for coverage. Here we will take the coordinates of TSS as example.

.. command-output:: gtftk select_by_key -i mini_real_noov_rnd_tx.gtf.gz -k feature -v transcript |  gtftk 5p_3p_coord > tss.bed
	:shell:

.. command-output:: gtftk mk_matrix -k 5 -u 5000 -d 5000 -i tss.bed -w 200 -l  H3K4me3,H3K79me,H3K36me3 ENCFF742FDS_H3K4me3_K562_sub.bw ENCFF947DVY_H3K79me2_K562_sub.bw ENCFF431HAA_H3K36me3_K562_sub.bw -o mini_real_single_nuc -c hg38.genome -t single_nuc
	:shell:


------------------------------------------------------------------------------------------------------------------


profile
~~~~~~~


Description: This command is used to create profil diagrams from a *mk_matrix* output. The two important arguments for
this command are *---group-by*, that defines the variable controling the set of colored lines and *---facet-var* that defines the variable controling the way the plot is facetted . Both *---group-by* and *---facet-var* should be set to one of *bwig*, *tx_classes* or *chrom*.


**Basic profiles**

A simple overlayed profile of all epigenetic marks around promoter. Here *---group-by* is, by default set to *bwig* and *---facet-var* is set to None. Thus a single plot with several lines corresponding to bwig coverage is obtained.


.. command-output:: gtftk profile -D -i mini_real_promoter.zip -o profile_prom -pf png -if example_01.png
	:shell:

.. image:: example_01.png
	:width: 75%


The same diagram is obtained if a bed file pointing to TSS was provided to *mk_matrix* and used in *single_nuc* mode.


.. command-output:: gtftk profile -i mini_real_single_nuc.zip -o profile_prom -pf png -if example_01a.png
	:shell:


.. image:: example_01a.png
	:width: 75%

Changing colors and applying color order can be done using the following syntax:


.. command-output:: gtftk profile -D -i mini_real_promoter.zip -c 'red,blue,violet' -d H3K79me,H3K4me3,H3K36me3 -o profile_prom -pf png -if example_01b.png
	:shell:


.. image:: example_01b.png
	:width: 75%


A subset of bigwig assessed for coverage can be selected for plotting. This is achieved using the *--subset-bwig* argument:

.. command-output:: gtftk profile -f bwig -g tx_classes -D -i mini_real_tx.zip  -fo  -o profile_tx -pf png -if example_01c.png  -fo -c 'red' -V 2  -w  -tl -e -lw 0.5 -u H3K4me3


.. image:: example_01c.png
	:width: 75%


Transcript coverage is obtained using the *mini_real_tx.zip* matrix. This provides a simple overlayed profile of all epigenetic marks along the transcript body extended in 5' and 3' regions:

.. command-output:: gtftk profile -D -i mini_real_tx.zip -o profile_tx -pf png -if example_02.png
	:shell:


.. image:: example_02.png
	:width: 75%


Almost the same but increasing the bins dedicated to upstream and dowstream regions (see *---bin-around-frac argument of *mk_matrix*).


.. command-output:: gtftk profile -D -i mini_real_tx_2.zip -o profile_tx -pf png -if example_03.png
	:shell:


.. image:: example_03.png
	:width: 75%


Note that the same is obtained when using user-defined regions (*i.e* when providing a bed as input corresponding to transcript coordinates).

.. command-output:: gtftk profile -D -i mini_real_user_def.zip -o profile_udef_4  -pf png -if example_04.png
	:shell:


.. image:: example_04.png
	:width: 75%


The same dataset used for plotting but adding a normalization step (*ranging*). When using *ranging* normalization, values are expressed as a percentage of the range between max and min value.


.. command-output:: gtftk profile -D -nm ranging -i mini_real_user_def.zip -o profile_udef_5  -pf png -if example_04b.png
	:shell:


.. image:: example_04b.png
	:width: 75%


Two examples using statistic 'max' and 2 differents values of '--upper-limit'.

.. command-output::  gtftk profile -D -i mini_real_promoter.zip -o profile_prom -pf png -if example_04_max_a.png  -V 2 -lw 1 -at 5 -s max -ul 1

.. image:: example_04_max_a.png
	:width: 75%



.. command-output::  gtftk profile -D -i mini_real_promoter.zip -o profile_prom -pf png -if example_04_max_b.png  -V 2 -lw 1 -at 5 -s max -ul 0.99

.. image:: example_04_max_b.png
	:width: 75%



**Faceted profiles**

Faceted plot of epigenetic profiles. The groups (i.e colors/lines) can be set to bwig classes and the facets to transcript classes. Things can be simply done by providing an additional file containing the transcript and their associated classes.


**Example:**


.. command-output:: gtftk profile -D -i mini_real_promoter.zip -f tx_classes -g bwig -fo -t tx_classes.txt -o profile_prom  -pf png -if example_05.png -e -V 2 -fc 2
	:shell:


.. image:: example_05.png
	:width: 100%


Alternatively, the groups can be set to chromosomes or transcript classes:


.. command-output:: gtftk profile -D -i mini_real_promoter.zip -g tx_classes -f bwig -fo -t tx_classes.txt -o profile_prom  -pf png -if example_06.png -V 2 -nm ranging
	:shell:


.. image:: example_06.png
	:width: 100%


.. command-output:: gtftk profile -D -i mini_real_promoter.zip -g chrom -f bwig -fo -t tx_classes.txt -o profile_prom  -pf png -if example_06b.png -V 2 -nm ranging
	:shell:


.. image:: example_06b.png
	:width: 100%


Note that facets may also be associated to epigenetic marks. In this case each the --group-by can be set to *tx_classes* or *chrom*.



.. command-output:: gtftk profile -D -i mini_real_tx_2.zip -g tx_classes -t tx_classes.txt -f bwig  -o profile_tx -pf png -if example_07.png  -fo -w -nm ranging
	:shell:


.. image:: example_07.png
	:width: 100%


.. command-output:: gtftk profile -D -i mini_real_tx_2.zip -g chrom -f bwig  -o profile_tx -pf png -if example_08.png  -fo -w -nm ranging
	:shell:


.. image:: example_08.png
	:width: 100%



**Theming**

The --theme argument controls plotnine theming.

.. command-output:: gtftk profile -th classic -D -i mini_real_promoter.zip -g bwig -f chrom  -o profile_prom  -c "#66C2A5,#FC8D62,#8DA0CB,#6734AF" -pf png -if example_09b.png
	:shell:

.. image:: example_09b.png
	:width: 100%


.. command-output:: gtftk profile -th seaborn -D -i mini_real_promoter.zip -g bwig -f chrom  -o profile_prom   -c "#66C2A5,#FC8D62,#8DA0CB,#6734AF" -pf png -if example_10.png
	:shell:

.. image:: example_10.png
	:width: 100%

.. command-output:: gtftk profile -th matplotlib -D -i mini_real_promoter.zip -g bwig -f chrom  -o profile_prom   -c "#66C2A5,#FC8D62,#8DA0CB,#6734AF" -pf png -if example_11.png
	:shell:

.. image:: example_11.png
	:width: 100%




**Playing with various commands**

It is also possible to use several of the previously seen commands to easily achieve more complexe analyses.
Here we will plot the epigenetic signal according to RNA-seq counts.


.. command-output:: gtftk join_attr -i mini_real_noov_rnd_tx.gtf  -j mini_real_counts_ENCFF630HEX.txt  -k gene_name -n counts | gtftk discretize_key -k  counts -n 6 -d count_levels -pu | gtftk tabulate -k transcript_id,count_levels -o tx_exp_classes.txt -Hun
	:shell:


.. command-output:: gtftk profile -D -i mini_real_tx.zip -o profile_tx -pf png -if example_12.png  -g tx_classes -t tx_exp_classes.txt -f bwig  -w -nm ranging -m viridis
	:shell:


.. image:: example_12.png
	:width: 100%


**Arguments:**

.. command-output:: gtftk profile -h
	:shell:


------------------------------------------------------------------------------------------------------------------


Commands from section 'miscellaneous'
---------------------------------------

col_from_tab
~~~~~~~~~~~~~~~~~~~~~~


**Description:** Select columns from a tabulated file based on their names.

**Example:**

.. command-output:: gtftk get_example | gtftk tabulate -k all |gtftk col_from_tab -c start,end,seqid | head -n 20
	:shell:

**Arguments:**

.. command-output:: gtftk col_from_tab -h
	:shell:


------------------------------------------------------------------------------------------------------------------


control_list
~~~~~~~~~~~~~~~~~~~~~~


**Description:** Returns a list of gene matched for expression based on reference values. Based on a reference gene list (or more generally IDs) this command tries to extract a set of other genes/IDs matched for signal/expression. The --reference-gene-file contains the list of reference IDs while the --inputfile contains a tuple gene/signal for all genes.

**Example:**

.. command-output:: gtftk control_list -i mini_real_counts_ENCFF630HEX.txt -r mini_real_control_1.txt -D -V 2 -s -l -p 1 -ju -if example_13.png -pf png
	:shell:


.. image:: example_13.png
	:width: 100%


**Arguments:**

.. command-output:: gtftk control_list -h
	:shell:


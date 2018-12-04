Commands from section 'Editing'
----------------------------------


add_prefix
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Add a prefix (or suffix) to one of the attribute value (*e.g.* gene_id)

**Example:**

.. command-output:: gtftk get_example| gtftk add_prefix -k transcript_id -t "novel_"| head -2
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

**Arguments:**

.. command-output:: gtftk del_attr -h
    :shell:


------------------------------------------------------------------------------------------------------------------

join_attr
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Add attributes from a file. This command can be used to import additional key/values into the gtf (e.g CPAT for coding potential, DESeq for differential analysis...). The imported file can be in 2 formats (2 columns or matrix):

- With a 2-columns file:

  - value for joining (transcript_id or gene_id...).
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

**Example:** With a matrix-like file.

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


**Example:** Add key/value to gene features.

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


**Arguments:**

.. command-output:: gtftk join_multi_file -h
    :shell:


------------------------------------------------------------------------------------------------------------------


discretize_key
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Create a new key by discretizing a numeric key. This can be helpful to create new classes of features on the fly.
The default is to create equally spaced interval. The intervals can also be created by computing the percentiles (-p) which will provide balanced classes most suitable generally.


**Example:** Let say we have the following matrix giving expression level of genes (rows) in samples (columns). We could join this information to the GTF and later choose to transform key *S1* into a new discretized key *S1_d*. We may apply particular labels to this factor using *-l*.


.. command-output:: gtftk get_example |  gtftk join_attr -j simple.join_mat -k gene_id -m | gtftk discretize_key -k S1 -d S1_d -n 2 -l A,B  | gtftk select_by_key -k feature -v gene
    :shell:

**Example:** We want to load RNA-seq data in the GTF (obtained through the get_example command) and discretize the expression values according to deciles (-p and -n set to 10). Classes will be labeled from A to J. The example below shows how balanced these classes will be.

.. seealso:: The *profile* command that could be used to asses the associated epigenetic marks of these 10 gene classes.


.. command-output:: gtftk get_example -d mini_real -f "*"
    :shell:

.. command-output:: gtftk get_example -d mini_real |  gtftk join_attr -H -j mini_real_counts_ENCFF630HEX.tsv -k gene_name -n exprs -t gene | gtftk discretize_key -k exprs -p -d exprs_class -n 10 -l A,B,C,D,E,F,G,H,I,J  | gtftk tabulate -k exprs_class -Hn | sort | uniq -c
    :shell:

**Arguments:**

.. command-output:: gtftk discretize_key -h
    :shell:


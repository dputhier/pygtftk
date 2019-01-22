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

get_5p_3p_coords
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Get the 5p or 3p coordinates for each feature (e.g TSS or TTS for a transcript).
Output is bed format.

**Example:** Get the 5p ends of transcripts and exons.

.. command-output:: gtftk get_example | gtftk get_5p_3p_coords -t transcript,exon -n transcript_id,gene_id,feature | head -n 5
	:shell:


**Arguments:**

.. command-output:: gtftk get_5p_3p_coords -h
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


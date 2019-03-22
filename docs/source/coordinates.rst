Commands from section 'coordinates'
-----------------------------------

In this section we will require the following datasets:

.. command-output:: gtftk get_example -q -d simple -f '*'
	:shell:


midpoints
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Get the genomic midpoint of each feature: genes, transcripts, exons or introns. Output is currently in bed format only.


**Example:** Get the midpoints of all transcripts and exons.

.. command-output:: gtftk midpoints -i simple.gtf -t transcript,exon -n transcript_id,feature | head -n 5
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

.. command-output:: gtftk get_5p_3p_coords  -i simple.gtf  -t transcript,exon -n transcript_id,gene_id,feature | head -n 5
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

.. command-output:: gtftk intergenic -i simple.gtf -c simple.chromInfo
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

.. command-output:: gtftk intronic -i simple.gtf | head -n 5
	:shell:


**Arguments:**

.. command-output:: gtftk intronic -h
	:shell:

------------------------------------------------------------------------------------------------------------------


splicing_site
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Compute the locations of donor and acceptor splice sites. This command will return a single position, which corresponds to the most 5' and/or the most 3' intronic region. If the gtf file does not contain exon numbering you can compute it using the
add_exon_nb command. The score column of the bed file contains the number of the closest exon relative to the splice site.

**Example:**

.. command-output:: gtftk add_exon_nb -i simple.gtf -k exon_nbr | gtftk splicing_site  -k exon_nbr| head
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

.. command-output:: gtftk shift -i simple.gtf  -s -10 -c simple.chromInfo | head -n 1
	:shell:


**Arguments:**

.. command-output:: gtftk shift -h
	:shell:


Commands from section 'annotation'
------------------------------------

In the example of this section we will need the following example files:

.. command-output:: gtftk get_example -q -d simple -f '*'
	:shell:


.. command-output:: gtftk get_example -q -d mini_real -f '*'
	:shell:


------------------------------------------------------------------------------------------------------------------


ologram
~~~~~~~~~~~~~~~~~~~~~~

**Description:** OLOGRAM -- OverLap Of Genomic Regions Analysis using Monte Carlo. Ologram annotates peaks
(in bed format) using (i) genomic features extracted from a GTF file (e.g promoter, tts, gene body, UTR...)
(ii) genomic regions tagged with particular keys/values in a GTF file (e.g. gene_biotype "protein_coding",
gene_biotype "LncRNA"...) or (iii) from a BED file (e.g. user-defined regions). Each couple peak file/region
is randomly shuffled across the genome (inter-region lengths are considered). Then the probability of intersection
under the null hypothesis (the peaks and this feature are independent) is deduced thanks to this Monte Carlo approach.
The program will return statistics for both the number of intersections and the total lengths (in basepairs) of all intersections.

.. warning:: The ologram examples below use 8 CPUs. Please adapt.



**Example:** Perform a basic annotation. We are searching whether H3K4me3 peaks tends to be enriched in some specific genomic elements.


.. command-output:: gtftk ologram -i mini_real.gtf.gz -p ENCFF112BHN_H3K4me3_K562_sub.bed -c hg38.genome -u 1500 -d 1500 -D  -if example_pa_01.pdf -k 8
	:shell:


.. raw:: html

  <br>
  <table>
  <tr>
  <td valign="top">
  <iframe src="_static/example_pa_01.pdf" title="your_title" align="top" width="500" height="620" width="50%" frameborder="0" scrolling="auto" target="Message">
  </iframe>
  </td>
  </tr>
  </table>
  <br>
  <br>

**Example:** Now we are using the gene_biotype key (note that a list of keys can be provided). This will tell us whether H3K4me3 tends to be located in particular transcripts (protein coding, LncRNAs...). The --no-basic-feature argument tells ologram not to test basic genomic elements (gene, transcripts...).

.. command-output:: gtftk ologram -i mini_real.gtf.gz -m gene_biotype -p ENCFF112BHN_H3K4me3_K562_sub.bed -c hg38.genome -D -n  -if example_pa_02.pdf -k 8
	:shell:


.. raw:: html

  <br>
  <table>
  <tr>
  <td valign="top">
  <iframe src="_static/example_pa_02.pdf" title="your_title" align="top" width="500" height="620" width="50%" frameborder="0" scrolling="auto" target="Message">
  </iframe>
  </td>
  </tr>
  </table>
  <br>
  <br>

**Example:** A more complex example where the key is created on the fly. Expression data are loaded as a novel key using the join_attr command and associated to gene features. This novel key (exprs) is then discretized to created 6 classes of genes with increasing expression (based on percentiles, -p) which are tested for enrichment in H3K4me3.

.. command-output:: gtftk join_attr -i mini_real.gtf.gz -H -j mini_real_counts_ENCFF630HEX.tsv -k gene_name -n exprs -t gene | gtftk discretize_key -k exprs -p -d exprs_class -n 6   | gtftk ologram -p ENCFF112BHN_H3K4me3_K562_sub.bed -c hg38.genome -D -n -m exprs_class -if example_pa_03.pdf -k 8
	:shell:


.. raw:: html

  <br>
  <table>
  <tr>
  <td valign="top">
  <iframe src="_static/example_pa_03.pdf" title="your_title" align="top" width="500" height="620" width="50%" frameborder="0" scrolling="auto" target="Message">
  </iframe>
  </td>
  </tr>
  </table>
  <br>
  <br>

**Example:** Using the add_exon_nb, we add the exon number transcript-wise (numbering from 5' to 3') and discretize this novel key into 5 classes tested for enrichment.

.. command-output:: gtftk add_exon_nb -k exon_nbr -i mini_real.gtf.gz | gtftk discretize_key -p -d exon_nbr_cat -n 5  -k exon_nbr | gtftk ologram -p ENCFF112BHN_H3K4me3_K562_sub.bed -c hg38.genome -D -n -m exon_nbr_cat -if example_pa_04.pdf -k 8
	:shell:


.. raw:: html

  <br>
  <table>
  <tr>
  <td valign="top">
  <iframe src="_static/example_pa_04.pdf" title="your_title" align="top" width="500" height="620" width="50%" frameborder="0" scrolling="auto" target="Message">
  </iframe>
  </td>
  </tr>
  </table>
  <br>
  <br>


**Arguments:**

.. command-output:: gtftk ologram -h
	:shell:



------------------------------------------------------------------------------------------------------------------


closest_genes
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Find the n closest genes for each transcript.

**Example:**

.. command-output:: gtftk closest_genes  -i simple.gtf -f
	:shell:


**Arguments:**

.. command-output:: gtftk closest_genes -h
	:shell:



------------------------------------------------------------------------------------------------------------------



overlapping
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Find transcripts whose body/TSS/TTS region extended in 5' and 3' (-u/-d) overlaps with any transcript from another gene. Strandness is not considered by default. Used --invert-match to find those that do not overlap. If --annotate-gtf is used, all lines of the input GTF file will be printed and a new key containing the list of overlapping transcripts will be added to the transcript features/lines (key will be 'overlapping_*' with * one of body/TSS/TTS). The --annotate-gtf and --invert-match arguments are mutually exclusive.


**Example:** Find transcript whose promoter overlap transcript from other genes.

.. command-output:: gtftk overlapping -i simple.gtf -c simple.chromInfo -t promoter -u 10 -d 10 -a    | gtftk select_by_key -k feature -v transcript | gtftk tabulate -k transcript_id,overlap_promoter_u0.01k_d0.01k | head
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
(Lepoivre C, BMC Genomics, 2013). The output is a GTF with an additional key
('divergent') whose value is set to '.' if the gene has no antisens transcript
in its promoter region. If the gene has an antisens transcript in its promoter
region the 'divergent' key is set to the identifier of the transcript whose tss
is the closest relative to the considered promoter. The tss to tss distance is
also provided as an additional key (dist_to_divergent).


**Example:** Flag divergent transcripts in the example dataset. Select them and produce a tabulated output.

.. command-output:: gtftk divergent -i simple.gtf -c simple.chromInfo -u 10 -d 10| gtftk select_by_key -k feature -v transcript | gtftk tabulate -k transcript_id,divergent,dist_to_divergent | head  -n 7
	:shell:

**Arguments:**

.. command-output:: gtftk divergent -h
	:shell:

------------------------------------------------------------------------------------------------------------------

convergent
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Find transcript with convergent tts. These transcripts will be defined here
as those whose tts region (defined by -u/-d) overlaps with the tts of
another gene in reverse/antisens orientation. The output is a GTF with an
additional key ('convergent') whose value is set to '.' if the gene has no
convergent transcript in its tts region. If the gene has an antisens transcript
in its tts region the 'convergent' key is set to the identifier of the
transcript whose tts is the closest relative to the considered tts.
The tts to tts distance is also provided as an additional key (dist_to_convergent).


**Example:** Flag divergent transcripts in the example dataset. Select them and produce a tabulated output.

.. command-output:: gtftk convergent -i simple.gtf -c simple.chromInfo -u 25 -d 25| gtftk select_by_key -k feature -v transcript | gtftk tabulate -k transcript_id,convergent,dist_to_convergent| head -n 4
	:shell:

**Arguments:**

.. command-output:: gtftk convergent -h
	:shell:

------------------------------------------------------------------------------------------------------------------

exon_sizes
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Add a new key to transcript features containing a comma-separated list of exon sizes.


**Example:**

.. command-output:: gtftk exon_sizes -i simple.gtf | gtftk select_by_key -t | gtftk tabulate -k transcript_id,exon_sizes
	:shell:

**Arguments:**

.. command-output:: gtftk exon_sizes -h
	:shell:

------------------------------------------------------------------------------------------------------------------


intron_sizes
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Add a new key to transcript features containing a comma-separated list of intron sizes.


**Example:**

.. command-output:: gtftk intron_sizes -i simple.gtf | gtftk select_by_key -t | gtftk tabulate -k transcript_id,intron_sizes
	:shell:

**Arguments:**

.. command-output:: gtftk intron_sizes -h
	:shell:


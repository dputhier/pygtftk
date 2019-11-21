Commands from section 'ologram'
------------------------------------


In the examples of this section, we will need the following example files:

.. command-output:: gtftk get_example -q -d simple -f '*'
	:shell:

.. command-output:: gtftk get_example -q -d mini_real -f '*'
	:shell:

.. command-output:: gtftk get_example -q -d hg38_chr1 -f '*'
	:shell:

.. command-output:: gtftk get_example -q -d ologram_1 -f '*'
	:shell:




**Note for contributors** : All files relevant to OLOGRAM are :
- *pygtftk/plugins/ologram.py* ; which is a wrapper.
- *pygtftk/plugins/merge_ologram_stats.py* ; a convenience function to merge OLOGRAM results.
- All files *in pygtftk/stats/intersect/* perform the calculations.



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


.. note:: The null hypothesis of the statistical test is:
	- H0: The regions of the query (--peak-file) are located independently of the reference (--inputfile or --more-bed) with respect to overlap.
	- H1: The regions of the query (--peak-file) tend to overlap the reference (--inputfile or --more-bed).


.. warning:: The ologram examples below use 8 CPUs. Please adapt.



**Example:** Perform a basic annotation. We are searching whether H3K4me3 peaks tends to be enriched in some specific genomic elements. The bars in
the bar plot diagram will be ordered according to 'summed_bp_overlaps_pvalue'.


.. command-output:: gtftk ologram -i hg38_chr1.gtf.gz -p ENCFF112BHN_H3K4me3_chr1.bed -c hg38_chr1.genome -u 1500 -d 1500 -D  -pf example_pa_01.pdf -k 8 -j summed_bp_overlaps_pvalue
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

<<<<<<< HEAD
=======


>>>>>>> 1.0.9
**Example:** Now we are using the gene_biotype key (note that a list of keys can be provided). This will tell us whether H3K4me3 tends to be located in particular transcripts (protein coding, LncRNAs...). The --no-basic-feature argument tells ologram not to test basic genomic elements (gene, transcripts...).

.. command-output:: gtftk select_by_key -i mini_real.gtf.gz -k gene_biotype -v protein_coding,lincRNA,antisense,processed_transcript  |  gtftk ologram  -m gene_biotype -p ENCFF112BHN_H3K4me3_K562_sub.bed -c hg38 -D -n  -pf example_pa_02.pdf -k 8 -j summed_bp_overlaps_pvalue
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


.. warning:: It may be important to consider the quality of the fit that is an indicator of the reliability of the p-value. This value is available in the tsv table produced by ologram. The fit quality may also be deplaced on the diagram using the -y/--display-fit-quality argument.


**Example:** A more complex example where the key is created on the fly. Expression data are loaded as a novel key using the join_attr command and associated to gene features. This novel key (exprs) is then discretized to created 6 classes of genes with increasing expression (based on percentiles, -p) which are tested for enrichment in H3K36me3.

.. command-output:: gtftk join_attr -i mini_real.gtf.gz -H -j mini_real_counts_ENCFF630HEX.tsv -k gene_name -n exprs -t exon | gtftk discretize_key -k exprs -p -d exprs_class -n 6  -u | gtftk ologram -p ENCFF119BYM_H3K36me3_K562_sub.bed -c hg38 -D -n -m exprs_class -pf example_pa_03.pdf -k 8 -j summed_bp_overlaps_pvalue
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

.. command-output:: gtftk add_exon_nb -k exon_nbr -i mini_real.gtf.gz | gtftk discretize_key -p -d exon_nbr_cat -n 5  -k exon_nbr | gtftk ologram -p ENCFF112BHN_H3K4me3_K562_sub.bed -c hg38 -D -n -m exon_nbr_cat -pf example_pa_04.pdf -k 8 -j summed_bp_overlaps_pvalue
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



merge_ologram_stats
~~~~~~~~~~~~~~~~~~~~~~

Several tsv files resulting from OLOGRAM analyses can be merged into a single diagram report using the merge_ologram_stats.

**Example:** For this example will will used the results obtained for 3 epigenetic marks on human chromosome 1.


.. command-output:: gtftk merge_ologram_stats H3K4me3_ologram_stats.tsv H3K79me2_ologram_stats.tsv H3K36me3_ologram_stats.tsv -l H3K4me3,H3K79me2,H3K36me3  -o merge_ologram_stats_01.pdf
	:shell:

.. raw:: html

  <br>
  <table>
  <tr>
  <td valign="top">
  <iframe src="_static/merge_ologram_stats_01.pdf" title="your_title" align="top" width="500" height="620" width="50%" frameborder="0" scrolling="auto" target="Message">
  </iframe>
  </td>
  </tr>
  </table>
  <br>



.. command-output:: gtftk merge_ologram_stats -h
	:shell:


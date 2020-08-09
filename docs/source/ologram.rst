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

.. command-output:: gtftk get_example -q -d simple_07 -f '*'
	:shell:

.. command-output:: gtftk get_example -q -d ologram_2 -f '*'
	:shell:



For more information about OLOGRAM and OLOGRAM-MODL, please see the appropriately titled papers in the Citing section.

More examples can be found in https://github.com/qferre/ologram_supp_mat and https://github.com/qferre/ologram-modl_supp_mat 
These contain example Snakemake workflows, that can be reused or from which commands can be extracted.

**Note for contributors** : To contribute to OLOGRAM, begin at *pygtftk/plugins/ologram.py* and unwrap function calls from there, to get a sense of how they interact. We have detailed comments to explain the role of every function.



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




**Example:** Perform a basic annotation. We are searching whether H3K4me3 peaks tends to be enriched in some specific genomic elements. The bars in the bar plot diagram will be ordered according to 'summed_bp_overlaps_pvalue'.


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


**Example:** We are now using the gene_biotype key (note that a list of keys can be provided). This will tell us whether H3K4me3 tends to be located in particular transcripts (protein coding, LncRNAs...). The --no-basic-feature argument tells ologram not to test basic genomic elements (gene, transcripts...).

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






**Example:** When not supplying a GTF, you can use --more-bed. The following example will look for pairwise enrichment of the file in input (p, here *query.bed* with the regions defined in --more-bed : here with *A.bed*, then with *B.bed*, then with *C.bed*.

.. command-output:: 
	gtftk ologram -ms 40 -mn 10 -p query.bed --more-bed A.bed B.bed C.bed -z -c hg38 -V 3 --force-chrom-peak --force-chrom-more-bed
  :shell:








ologram (multiple overlaps)
~~~~~~~~~~~~~~~~~~~~~~

While previously we computed paiwise enrichment (ie. Query+A, Query+B ...) , It is also possible to use the **OLOGRAM-MODL** Multiple Overlap Dictionary Learning) plugin to find multiple overlaps (ie. between n>=2 sets) enrichment (ie. Query+A+B, Query+A+C, ...) in order to highlight combinations of genomic regions, such as Transcriptional Regulator complexes. 

This is done only on custom regions supplied as BEDs supplied with the `--more-bed` argument. In most cases you may use the --no-gtf argument and only pass the regions of interest.


For statistical reasons, we recommend shuffling across a relevant subsection of the genome only (ie. enhancers only) using --bed-excl or --bed-incl to ensure the longer combinations have a reasonable chance of being randomly encountered in the shuffles.


**MODL itemset mining algorithm** :By default, OLOGRAM-MODL will compute the enrichment of all n-wise combinations that are encountered in the real data it was passed. This however can add up to 2**N combinations and make the result hard to read. As such, we also give the option to use a custom itemset mining algorithm on the true overlaps to identify interesting combinations. 

In broad strokes, this custom algorithm MODL (Multiple Overlap Dictionary Learning) will perform many matrix factorizations on the matrix of true overlaps to identify relevant correlation groups of genomic regions. Then a greedy algorithm based on how much these words improve the reconstruction will select the utmost best words. MODL is only used to filter the output of OLOGRAM : once it returns a list of interesting combination, OLOGRAM will compute their enrichment as usual, but for them only. Each combination is of the form [Query + A + B + C] where A, B and C are BED files given as --more-bed. You can also manually specify the combinations to be studied with the format defined in OLOGRAM notes (below).

Unlike classical association rules mining algorithms, this focuses on mining relevant bio complexes/clusters and correlation groups (item sets), and you should not request more than 20-30 combinations. This itemset mining algorithm is a work-in-progress. If you want the enrichment of all combinations, ignore it. To use MODL, use the --multiple-overlap-max-number-of-combinations argument.


**Exact combinations ** : By default, OLOGRAM will compute "inexact" combinations, meaning that when encountering an overlap of [Query + A + B + C] it will count towards [A + B + ...]. For exact intersections (ie. [Query + A + B + nothing else]), set the --multiple-overlap-target-combi-size flag to the number of --more-bed plus one. You will know if the combinations are computed as inexact by the "..." in their name in the result file. Intersections not including the query file are discarded.





**Example:**

.. command-output:: 

  gtftk ologram -z -c simple_07.chromInfo -p simple_07_peaks.bed       # The query (-p) is the file to compare against.
      --more-bed simple_07_peaks.1.bed simple_07_peaks.2.bed           # List of files to compare with
      # --more-bed `ls -d ./data/*`                                    # This should work instead if all your files are in the 'data' subdirectory
      -o results --force-chrom-peak --force-chrom-more-bed  
      -V 3 -k 8 -mn 10 -ms 10                                          # Verbosity, threads, number and size of minibatches
      --more-bed-multiple-overlap                                      # Use multiple overlaps on the --more-bed
      --multiple-overlap-max-number-of-combinations 10                 # OPTIONAL ARGUMENT. Use MODL to restrict to THIS MANY combinations (optional)
      --multiple-overlap-target-combi-size 3                           # OPTIONAL ARGUMENT. Combis restricted to this size. Also Explain exact (optional)
      --multiple-overlap-custom-combis test_combis.txt                 # OPTIONAL ARGUMENT. Will bypass the selection by the previous two arguments and work only on the combinations defined in this file.
  :shell:


.. raw:: html

  <br>
  <table>
  <tr>
  <td valign="top">
  <iframe src="_static/example_ologram_modl.pdf" title="your_title" align="top" width="500" height="620" width="50%" frameborder="0" scrolling="auto" target="Message">
  </iframe>
  </td>
  </tr>
  </table>
  <br>
  <br>





**MODL algorithm API:** MODL can also be used independantly as a combination mining algorithm. You data needs to be a matrix with one line per transaction and one column per element.

For more details, see code comments.

Here is an example:



``` python
from pygtftk.stats.intersect.dict_learning import Modl, test_data_for_modl
flags_matrix = test_data_for_modl(nflags = 1000, number_of_sets = 6, noise = 0.1, cor_groups = [(0,1),(0,1,2,3),(4,5)])
combi_miner = Modl(flags_matrix, 
       multiple_overlap_target_combi_size = -1,    # Limit the size of the combinations
       multiple_overlap_max_number_of_combinations = 3,    # How many words to find ?
       nb_threads = 1,
       step_1_factor_allowance = 2)    # How many words to ask for in each step 1 rebuilding
interesting_combis = combi_miner.find_interesting_combinations()   
```



For more details about usage and implementation, please read the notes below :

**Arguments:**

.. command-output:: gtftk ologram -h
	:shell:





ologram_merge_stats
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Several tsv files resulting from *OLOGRAM* analyses can be merged into a single diagram report using the merge_ologram_stats.

**Example:** For this example we will used the results obtained for 3 epigenetic marks on human chromosome 1.

.. command-output:: gtftk ologram_merge_stats H3K4me3_ologram_stats.tsv H3K36me3_ologram_stats.tsv H3K79me2_ologram_stats.tsv -o merge_ologram_stats_01.pdf --labels H3K4me3,H3K36me3,H3K79me2
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
  <br>

This also works with OLOGRAM-MODL results, since they follow the same basic format of one element/combination per line.

**Arguments:**

.. command-output:: gtftk ologram_merge_stats -h
	:shell:




ologram_modl_treeify
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Visualize n-wise enrichment results (OLOGRAM-MODL) as a tree of combinations. Works on the result (tsv file) of an OLOGRAM analysis called with --more-bed-multiple-overlap.

We recommend this representation. The tsv file can be edited before passing it to the command, for example by keeping only the combinations you are interested in.

On the graph, S designated the total number of basepairs in which this combinations is encountered in the real data. Fold change gives the ratio with the number of basepairs in the shuffles, with the associated Negative Binomial p-value.

.. command-output:: gtftk ologram_modl_treeify -i multiple_overlap_trivial_ologram_stats.tsv -o ./results/treeified.pdf -l ThisWasTheNameOfTheQuery
	:shell:

.. raw:: html

  <br>
  <table>
  <tr>
  <td valign="top">
  <iframe src="_static/example_ologram_treeify.pdf" title="your_title" align="top" width="500" height="620" width="50%" frameborder="0" scrolling="auto" target="Message">
  </iframe>
  </td>
  </tr>
  </table>
  <br>
  <br>

.. command-output:: gtftk ologram_modl_treeify -h
	:shell:




ologram_merge_runs
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Merge several runs of OLOGRAM into a single run, by treating each a "superbatch" of shuffles.

OLOGRAM remembers all intersections occuring inside all minibatches, so as to calculate statistics. If you are using a large number of shuffles and/or very large files, this may cost a lot of RAM. In practice, you will seldom need more than 100 shuffles. But optionally, if you require increased precision, you can run OLOGRAM several times, treat each run as a "batch of batches" and merge and recalculate stats on the merged superbatch automatically using this command.

```bash
# Make several OLOGRAM runs
N_RUNS = 100
for i in {1..$N_RUNS}
do
   ologram ...
done

# Merge those runs
gtftk ologram_merge_runs --inputfiles `ls ./results/*.tsv` -o ./merged_batches_result.tsv -V 3
```

Other commands such as ologram_modl_treeify can now be called on the resulting tsv, which respects the OLOGRAM format.


.. command-output:: gtftk ologram_merge_runs -h
	:shell:


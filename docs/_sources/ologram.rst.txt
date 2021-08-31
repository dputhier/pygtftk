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

More examples can be found in <https://github.com/qferre/ologram_supp_mat> and more recently at <https://github.com/qferre/ologram-modl_supp_mat> 
These contain example Snakemake workflows, that can be reused or from which commands can be extracted.

Most of the commands presented in this section are demonstrated in the *ologram-modl_supp_mat* Git, along with certain perspectives.

**Note for contributors** : To contribute to OLOGRAM, begin at *pygtftk/plugins/ologram.py* and unwrap function calls from there, to get a sense of how they interact. We have detailed comments to explain the role of every function. *A detailed table with the role of each file is presented at the end of this document.*



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


.. warning:: The ologram examples below use 8 CPUs. Please adapt the number of threads.




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






**Example:** When not supplying a GTF, you can use --more-bed. The following example will look for pairwise enrichment of the file in input (p, here *query.bed* with the regions defined in --more-bed : here query with *A.bed*, then query with *B.bed*, then query with *C.bed*.

.. code-block:: bash

	gtftk ologram -ms 40 -mn 10 -p query.bed --more-bed A.bed B.bed C.bed -z -c hg38 -V 3 --force-chrom-peak --force-chrom-more-bed









ologram (multiple overlaps)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is also possible to use the **OLOGRAM-MODL** Multiple Overlap Dictionary Learning) plugin to find multiple overlaps (ie. between n>=2 sets) enrichment (ie. Query+A+B, Query+A+C, ...) in order to highlight combinations of genomic regions, such as Transcriptional Regulator complexes. 

.. note:: The null hypothesis of the statistical test is:
	- H0: Considering the genomic regions in the query set (--peak-file) and in the reference sets (--more-bed), the regions in one set are located independently of the regions in any another set. They are not assumed to be uniformly distributed, we keep inter-region lengths.
              
This is done only on custom regions supplied as BEDs supplied with the `--more-bed` argument. In most cases you may use the --no-gtf argument and only pass the regions of interest.

For statistical reasons, we recommend shuffling across a relevant subsection of the genome only (ie. enhancers only) using --bed-excl or --bed-incl. This ensures the longer combinations have a reasonable chance of being randomly encountered in the shuffles. Conversely, if you do not filter the combinations, keep in mind that the longer ones may be enriched even though they are present only on a few base pairs, because at random they would be even rarer. As such, we recommend focusing comparisons on combinations of similar order (number of sets).

**Exact combinations:** By default, OLOGRAM will compute "inexact" combinations, meaning that when encountering an overlap of [Query + A + B + C] it will still count as an observation of [Query + A + B + ...] (meaning "Query + A + B + any other set"). For exact intersections (ie. [Query + A + B + nothing else]), set the --exact flag to True. You will know if the combinations are computed as inexact by the '...' in their name in the result file. 

In any case, only intersections with the query are counted. ie. Query+A+B is counted, but A+B+C is not.

With inexact combinations, if A+B is very enriched and C is depleted, A+B+C will be enriched. It is more interesting to look at C's contribution to the enrichment. Relatedly, longer combinations are usually more enriched since they involve more theoretically independant sets. Relatedly, you should compare the enrichments of combinations of similar orders (number of sets in the combinations) since longer combinations tend to be more enriched under (H_0).



**Simple example:**

Comparing the query (-p) against two other BED files, analyzing multiple overlaps.

.. command-output:: gtftk ologram -z -w -q -c simple_07.chromInfo -p simple_07_peaks.bed --more-bed simple_07_peaks.1.bed simple_07_peaks.2.bed --more-bed-multiple-overlap
  :shell:


**Detailed example:**

.. code-block:: bash

  gtftk ologram -z -c simple_07.chromInfo -p simple_07_peaks.bed     # The query (-p) is the file to compare against.
    --more-bed simple_07_peaks.1.bed simple_07_peaks.2.bed           # List of BED files giving the region sets to compare with. TIP: You can use  --more-bed `ls -d ./data/*` if all your files are in the 'data' subdirectory
    -o results --force-chrom-peak --force-chrom-more-bed  
    -V 3 -k 8 -mn 10 -ms 10                                          # Verbosity, threads, number and size of minibatches
    --more-bed-multiple-overlap                                      # Toggle the computation of multiple overlaps on the --more-bed
    --exact                                                          # OPTIONAL ARGUMENT. If present, an observation of A+B+C will not count as an observation of A+B.
    --multiple-overlap-max-number-of-combinations 10                 # OPTIONAL ARGUMENT. Use MODL to restrict to this many combinations.
    --multiple-overlap-target-combi-size 3                           # OPTIONAL ARGUMENT. Combis mined longer than this size will not be shown.
    --multiple-overlap-custom-combis test_combis.txt                 # OPTIONAL ARGUMENT. Will bypass the selection by the previous two arguments and work only on the combinations defined in this file.
    --keep-intact-in-shuffling 0,1                                   # BETA - OPTIONAL ARGUMENT. Gives the positions of the files in --more-bed that will be kept fixed in shuffling.

See the result of `gtftk ologram -h` below for more detailed informations about the arguments' formats.


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


As the computation of multiple overlaps can be RAM-intensive, if you have a very large amount of candidate genomic feature sets (hundreds) we recommend selecting less candidates among them first by running a pairwise analysis.


**MODL itemset mining algorithm:** By default, OLOGRAM-MODL will compute the enrichment of all n-wise combinations that are encountered in the real data it was passed. This however can add up to 2**N combinations and make the result hard to read. Furthermore, in biological data noise is a real problem and can obscure the relevant combinations. As such, we also give the option to use a custom itemset mining algorithm on the true overlaps to identify interesting combinations. Another possibility is to instead manually pass a text file containg the combinations you want to study



Itemset mining details
======================

In broad strokes, the custom itemset algorithm MODL (Multiple Overlap Dictionary Learning) will perform many matrix factorizations on the matrix of true overlaps to identify relevant correlation groups of genomic regions. Then a greedy algorithm based on how much these words improve the reconstruction will select the utmost best words. MODL is only used to filter the output of OLOGRAM : once it returns a list of interesting combination, OLOGRAM will compute their enrichment as usual, but for them only. Each combination is of the form [Query + A + B + C] where A, B and C are BED files given as --more-bed. You can also manually specify the combinations to be studied with the format defined in OLOGRAM notes (below).

Unlike classical association rules mining algorithms, this focuses on mining relevant biological complexes/clusters and correlation groups (item sets). As such, we do not recommend asking for more than 20-50 combinations to keep the running time reasonable and keep the found combinations still relevant.

As a matrix factorization based algorithm, it is designed to be resistant to noise which is a known problem in biological data. Its goal is to extract meaningful frequent combinations from noisy data. As a result however, it is biased in favor of the most abundant combinations in the data, and may return correlation groups if you ask for too few words (ie. if AB, BC and AC are complexes, ABC might be returned).

This itemset mining algorithm is a work-in-progress, and optional . Whether you use MODL will not change the results for each combination, it only changes which combinations are displayed. If you want the enrichment of all combinations, ignore it. To use MODL, use the --multiple-overlap-max-number-of-combinations argument.

MODL is mostly needed when the list of -\-more-bed is very long and you do not want to filter the results manually, and when you are working with noisy data which could obfuscate the interesting combinations. It is also possible to bypass it and provide a custom list of combinations to be considered.

 

**MODL algorithm API:** MODL can also be used independantly as a combination mining algorithm. 

This can work on any type of data, biological or not, that respects the conventional formatting for lists of transactions: the data needs to be a matrix with one line per transaction and one column per element. For example, if you have three possible elements A, B and C, a line of [1,0,1] means a transaction containing A and C.

For a factor allowance of k and n final queried words, the matrix will be rebuilt with k*n words in step 1. MODL will discard combinations rarer than 1/10000 occurences to reduce computing times. It will also reduce the abundance of all unique lines in the matrix to their square roots to reduce the emphasis on the most frequent elements. However, the latter can magnify the impact of the noise as well and can be disabled when using the manual API. To de-emphasize longer words, which can help in this case, we normalize words by their summed square in step 2.

If you are passing a custom error function, it must have this signature: `error_function(X_true, X_rebuilt, encoded, dictionary)`. X_true is the real data, and X_rebuilt is the reconstruction to evaluate.
encoded is the encoded version (U) which in our case is used to assess sparsity, while dictionary (V) is the matrix with one atom of the dictionaty per row (not used by default). Note that the dictionary is passed before MODL performs any normalization on it.  All are NumPy matrices.


.. note:: An example of custom loss we recommend is: selecting the combinations (of reference sets) that best predict the query set using a Naive Bayes classifier. This is not yet implemented, but a fully functional example is available at <https://github.com/qferre/ologram-modl_supp_mat/blob/master/scripts/modl_perspective.py> as a Python script. To use it, simply replace the filepaths at the beginning with the paths to your own files, and run the script. The order in the selection will be the same as the order you gave in the script, not alphabetical. You can then run OLOGRAM without MODL, or pass the custom selection you just computed.




**For more details, see code comments.**

Here is an example:

.. code-block:: python

  from pygtftk.stats.intersect.modl.dict_learning import Modl, test_data_for_modl
  flags_matrix = test_data_for_modl(nflags = 1000, number_of_sets = 6, noise = 0.1, cor_groups = [(0,1),(0,1,2,3),(4,5)])

  from pygtftk import utils
  utils.VERBOSITY = 2 # Ensure DEBUG messages are shown


  # Simple example of custom error function
  def my_error_function (X_true, X_rebuilt, encoded, dictionary):  return np.sum(X_rebuilt - X_true)

  combi_miner = Modl(flags_matrix, 
    multiple_overlap_target_combi_size = -1,            # Limit the size of the combinations
    multiple_overlap_max_number_of_combinations = 3,    # How many words to find ?
    nb_threads = 1,
    step_1_factor_allowance = 2,                        # (Defaults to 2) How many words to ask for in each step 1 rebuilding, as a multiplier of multiple_overlap_max_number_of_combinations.
    error_function = None,                              # (OPTIONAL) Custom error function in step 2
    smother = True,                                     # (Defaults to True) Should the smothering (quadratic reduction of abundance) be applied ?
    normalize_words = True,                             # (Defaults to True) Normalize words by their summed squared in step 2 ?
    step_2_alpha = None,                                # (OPTIONAL) Override the alpha (sparsity control) used in step 2.
    discretization_threshold = 0                        # (Defaults to 1) Discretization threshold D : in each atom, elements below D*maximum_for_this_atom will be discarded.
    step_1_alphas = None)                               # (OPTIONAL) Override the list of alphas used in step 1 (should be a list)
  interesting_combis = combi_miner.find_interesting_combinations()   


For more details about usage and implementation, please read the notes below.

**Arguments:**

.. command-output:: gtftk ologram -h
	:shell:



**Manual intersection computing:** To manually compute an overlap matrix between any number of BED files, the following Python code can be used.

.. code-block:: python

  import pybedtools
  import numpy as np
  from pygtftk.stats.intersect.overlap_stats_compute import compute_true_intersection

  # Some example paths
  QUERY_PATH = "./input/query.bed"
  MORE_BED_PATHS = ["./input/A.bed", "./input/B.bed", "./input/C.bed"]
  EXCL_PATH = "./exclusion.bed"

  # Register the BED files as pybedtools.BedTool objects
  bedA = pybedtools.BedTool(QUERY_PATH)
  bedsB = [pybedtools.BedTool(bedfilepath).sort().merge() for bedfilepath in MORE_BED_PATHS] # Sort and merge for the bedsB

  # OPTIONAL - Exclude some regions from the BEDs
  bed_excl = pybedtools.BedTool(EXCL_PATH)
  bedA = read_bed.exclude_concatenate(bedA, bed_excl)
  bedsB = [read_bed.exclude_concatenate(bedB, bed_excl) for bedB in bedsB]

  # Use our custom intersection computing algorithm to get the matrix of overlaps
  true_intersection = compute_true_intersection(bedA, bedsB)
  flags_matrix = np.array([i[3] for i in true_intersection])

  # If desired, run MODL or any other algorithm on this
  my_algorithm.process(flags_matrix)
  # See code block above for a MODL example

The resulting flags_matrix is a NumPy array that can be edited, and on which MODL can be run. It is also possible to run any itemset miner you wish on this matrix. An implementation of apriori is provided in the `pygtftk.stats.intersect.modl.apriori.Apriori` class.

Note that by definition, in this intersections' matrix only regions where at least two sets are open are given. Regions where a single set was open will not be present.
If you want a matrix of all contiguous elements where at least one set is open, and not just intersections, you may opt to instead use as "query" (bedA) a BED file covering all the chromosomes in the genome (e.g. if your genome has only 2 chromosomes of length 100 each, this "query" file would be "chr1 0 100 \n chr2 0 100"). 
To have predictable binning based on length in the final matrix instead of one line per intersection, you may also subdivide fake "query "chr1 0 100" region into bins of, say, 10 bp instead: "chr1 0 10 \n chr1 11 20\n ...".

Since the results of MODL only depend on the true intersections and not on the shuffles, you can run also MODL with 1 shuffle or on a manually computed matrix as above to pre-select interesting combinations, and then run the full analysis on many shuffles. We then recommend selecting the combinations that interest you in the resulting tsv file, using MODL's selection as a starting point and adding or removing some combinations based on your own needs (eg. adding all the highest fold changes, or all particular combinations containing the Transcription Factor X that you are studying).



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

Cases without a p-value diamond mean it was NaN. It usually means was too rare to be encountered in the shuffles.

An example of use case for this tool would be to compare between different cell lines, or to slop (extend) your query regions by different lengths and compare the enrichment to find at which distance of each other several sets are on average.

**Arguments:**

.. command-output:: gtftk ologram_merge_stats -h
	:shell:



ologram_modl_treeify
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Visualize n-wise enrichment results (OLOGRAM-MODL) as a tree of combinations. Works on the result (tsv file) of an OLOGRAM analysis called with --more-bed-multiple-overlap. On the graph, S designated the total number of basepairs in which this combinations is encountered in the real data. Fold change gives the ratio with the number of basepairs in the shuffles, with the associated Negative Binomial p-value.

This recommended representation is useful to find master regulators, by showing which additions to a combinations increase its enrichment, and allowing to see whether overlaps that contain the element X also contain the element Y (looking at how a child combination accounts for the S of its parent in an inexact counting).

P-values of NaN (-1 in the original tsv) are due to poor fitting. They are mostly present in high order combinations, that were so rare that they are not encountered in the shuffles even once. We also recommend discarding the rarest combinations found on such a very small number of basepairs that they are unlikely to be biologically significant. This is mostly relevant when you have many sets (k >= 5) since longer combinations will often be enriched through sheer unlikelihood. To that effect, there is a parameter to display only the combinations with the highest S.

The tsv result file can be edited before passing it to the command, for example by keeping only the combinations you are interested in. 
You can either (1) run OLOGRAM-MODl with no filtering and get a tree of all combinations, (2) use MODL to get a pre-selection that can be tailored, or (3) take the run with all combinations from the possibility 1 and use the -t argument to take the most frequent combinations.

.. command-output:: gtftk ologram_modl_treeify -i multiple_overlap_trivial_ologram_stats.tsv -o treeified.pdf -l ThisWasTheNameOfTheQuery
	:shell:

.. raw:: html

  <br>
  <table>
  <tr>
  <td valign="top">
  <iframe src="_static/treeified.pdf" title="your_title" align="top" width="500" height="620" width="50%" frameborder="0" scrolling="auto" target="Message">
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

OLOGRAM remembers all intersections occuring inside all minibatches, so as to calculate statistics. If you are using a large number of shuffles and/or very large files, this may cost a lot of RAM. In practice, you will seldom need more than 100-200 shuffles. But optionally, if you require increased precision, you can run OLOGRAM several times, treat each run as a "batch of batches" and merge and recalculate stats on the merged superbatch automatically using this command.

Around 100-200 shuffles is usually enough to robustly fit a Negative Binomial distribution. In terms of precision, a Negative Binomial mean under 1/100 (meaning this combination was not seen at least once in 100 shuffles) would not mean much anyways. 

.. code-block:: bash

  # Make several OLOGRAM runs
  N_RUNS = 100
  for i in {1..$N_RUNS}
  do
    gtftk ologram ...     # Replacing this with a complete OLOGRAM command
  done

  # Merge those runs
  gtftk ologram_merge_runs --inputfiles `ls ./results/*.tsv` -o ./merged_batches_result.tsv -V 3


Other commands such as ologram_modl_treeify can now be called on the resulting tsv, which respects the OLOGRAM format.

.. command-output:: gtftk ologram_merge_runs -h
	:shell:






Notes
~~~~~~~~~~~~~~~~~~~~~~

*This section contains more specific notes about the use and interpetation of OLOGRAM*.

-- The goal of the minibatches is to save RAM. You should increase the number of minibatches, instead of their size.

-- If -\-more-keys is used additional region sets will be tested based on the associated key value. As an example, if -\-more-keys is set to the 'gene_biotype' (a key generally found in ensembl GTF), the region related to 'protein_coding', 'lncRNA' or any other values for that key will be retrieved merged and tested for enrichment.

-- For statistical reality reasons, with multiple sets the expected overlaps for the longer combinations (A+B+C+D+... when they are all independant) can be very low. As a result, longer combinations tend to be more enriched: this should be kept in mind when comparing enrichment values between combinations of a different order. 

This is especially true when the total genomic coverage of the sets is low. We recommend instead shuffling only across a biologically relevant subsection of the genome (with -\-bed-incl) : for example, if studying  Transcriptional Regulators, shuffling only on inferred Cis Regulatory Modules or promoters.
                          If the shuffling is restricted to a sub-genome, and features outside are discarded. In essence it mostly means switching to a smaller genome. Of course, since the shuffling is done only here, (H_0) becomes ‘... the features are independent and can only be located on the sub-genome’. This bears mentioning. In practice, this means shuffling only across shortened chromosomes.

Our Negative Binomial model helps alleviate this problem. Even so, if a combination is so rare that it is not encoutered even once in the shuffles, it will have a p-value of NaN. Furthermore, if C is depleted with query but always present with A and B, and A and B are enriched themselves, A+B+C will be enriched.

-- BETA - When using --more-bed (and only that), you can give a list of bed files that should be kept fixated during the shuffles using the --keep-intact-in-shuffling argument.

-- RAM will be the biggest limiting factor. While 100 total shuffles should be enough to fit a Negative Binomial distribution in most cases, if needed try running more batches of fewer shuffles instead of the other way around. The alternative is running them independantly and merging them afterwards with *ologram_merge_runs*.

-- If you have many (30+) BED files in --more-bed, consider running a pairwise analysis first to divide them in groups of 10-20, and study the multiple overlaps within those groups. This is also more likely to be biologically significant, as for example Transcription Factor complexes usually have 2-8 members.

-- We recommend running the ologram_modl_treeify plugin on the resulting tsv file if you use multiple overlaps.


-- Our Negative Binomial model is only an approximation for the underlying true distribution, which is likely close to a Beta Binomial. For instance, the Neg. Binom. approximation fails with too few regions in the sets (at least 1K), and will likely slightly overestimate the p-values in other cases. However, precision is usually good for even very significant p-values, dropping only at the very significant level (<1E-5), hence there is only a very small risk of false positives. Also, even if they are overestimated, the order of p-values is unchanged (as a Neg. Binom. is a special case of Beta) meaning if a combination 1 has a higher Neg. Binom. p-value than combination 2, its true p-value is also likely higher than the p-value of 2.

The Neg. Binom. is still the better option, as fitting the proper distribution (approximated as Beta) is more difficult. As such, an ad-hoc p-value based on the Beta distribution is given, but it will only better than the Neg. Binom. on massive numbers of shuffles (thousands). We also added the empirical p-value as a new column (ie. number of shuffles in which a value as extreme is observed) if you believe the model to be inadequate.

-- Our model rests upon certain assumptions (ie. exchangeable variables, sufficient nb. of regions, etc.). The null hypothesis can be rejected if any assumption is rejected, or merely because the approximation holds only asymptotically. The fitting test is the key for that: if, when performing the shuffles, it is found that the distribution of S under our shuffling model does not follow a Neg. Binom., it will be said. Then if the hypothesis is rejected (low p-val) but the fitting was good, it is then reasnobale to assume the combination is enriched. Admittedly, the fitting test does not test the tails of the distribution, but it shows if the general shape is close enough.


------------------------------------------------------------------------------------------------------------------


OLOGRAM file structure
~~~~~~~~~~~~~~~~~~~~~~

Below is a detailed list of the source code files of OLOGRAM-MODL, with their roles. All paths are relative to the root of the *pygtftk* Git repository. The "Plugin" group designates plugins that can be called directly from the command line. A file extension of "pyx" designates a Cython file, "py" a Python file, and "cpp" a C++ file.


.. list-table:: OLOGRAM-MODL files.
   :widths: 10 40 50
   :header-rows: 1

   * - Group
     - File path
     - Description
   * - Plugin
     - pygtftk/plugins/ologram.py
     - *Root file.* Parses the arguments and calls the other functions.
   * - Utility
     - docs/source/ologram.rst
     - Documentation source.
   * - Root
     - pygtftk/stats/intersect/overlap_stats_shuffling.py
     - *Main function*. Called directly by the *ologram.py* plugin. All other functions calls are descended from this one.
   * - Root
     - pygtftk/stats/intersect/overlap_stats_compute.py
     - Functions to compute overlap statistics on (shuffled) region sets. 
   * - Algorithm
     - pygtftk/stats/intersect/create_shuffles.pyx
     - Shuffle BED files and generate new "fake" BED files.
   * - Algorithm
     - pygtftk/stats/intersect/overlap/overlap_regions.pyx
     - Compute the overlaps between two sets of genomic regions, supporting multiple overlaps.
   * - Utility 
     - Turn a BED file into a list of intervals.
     - pygtftk/stats/intersect/read_bed/read_bed_as_list.pyx
   * - Utility
     - pygtftk/stats/intersect/read_bed/exclude.cpp
     - Exclude certain regions from a set to create concatenated sub-chromosomes.
   * - Utility
     - pygtftk/stats/multiprocessing/multiproc.pyx
     - Helper functions and structures for multiprocessing.
   * - Statistics
     - pygtftk/stats/negbin_fit.py
     - Utility functions relative to the negative binomial distribution, including verifying its good fit.
   * - MODL
     - pygtftk/stats/intersect/modl/dict_learning.py
     - Contains the MODL algorithm, an itemset mining algorithm described in this paper.
   * - MODL
     - pygtftk/stats/intersect/modl/subroutines.py
     - Subroutines of the MODL algorithm. Those are pure functions and can be used independently.
   * - Utility
     - pygtftk/stats/intersect/modl/tree.py
     - A graph-based representation of combinations of elements.
   * - Plugin
     - pygtftk/plugins/ologram_merge_runs.py
     - Merge a set of OLOGRAM runs into a single run and recalculates statistics based on it.
   * - Plugin
     - pygtftk/plugins/ologram_merge_stats.py
     - Merge a set of OLOGRAM outputs calculated on different queries into a single output, preserving labels. Build a heatmap from the results.
   * - Plugin
     - pygtftk/plugins/ologram_modl_treeify.py
     - Turns a result of OLOGRAM-MODL multiple overlap (tsv file) in a tree for easier visualisation.
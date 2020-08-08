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



Need those too no ?

.. command-output:: gtftk get_example -q -d simple_07 -f '*'
	:shell:

.. command-output:: gtftk get_example -q -d ologram_2 -f '*'
	:shell:





For more information about OLOGRAM and OLOGRAM-MODL, please see the following papers :
	- paper of ologram, Ferré et al., 2019, add link once published
	- paper of ologram-modl, Ferré et al., 2020?, add link once submitted and published




VERY IMPORTANT :
	More examples can be found in https://github.com/qferre/ologram_supp_mat as a Snakemake workflow, or you can extract commands
  Also add the github of ologoram modl





**Note for contributors** : All files relevant to OLOGRAM are :
- *pygtftk/plugins/ologram.py* ; which is a wrapper.
- *pygtftk/plugins/ologram_merge_stats.py* ; a convenience function to merge OLOGRAM results.
- All files *in pygtftk/stats/intersect/* perform the calculations.
REMOVE THIS ?

Start at ologram.py and unwrap function calls from there, to get a sense of how they interact.
Main one is overlap_stats_shuffling.compute_overlap_stats, which itself calls the functions in overlap_stats_compute, which then calls functions from other modules.

We have detailed comments to explain the role of every function





Add links to the paper for the details ! Note to do so in my Zim ! So I do not have to re-explain all the algos here.



Add an example with ologram_merge_stats

Add an exemple of ologram with multiple overlap with custom combis, and one with dict learning. Add a short explanation, but mostly link to the paper.
If there is algorithmical novelty in my dict learning approach, add details










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



There is  MULTIPROCESS BY BATCH ! EXPLAIN IT TO THEM ! Mostly that it will also increase RAM usage, predictably.








Examples of use are available at : <give the link to the two githubs for ologram supp mat et ologram modl supp mat>







More details, see the NOTES !!! They should be below with the argument












Add direct link to ologram.py source on github

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


























**Example:** When not supplying a gtf, using --more-bed

.. command-output:: 
	:shell:

	gtftk ologram -ms 40 -mn 10 -p query.bed \
            --more-bed A.bed B.bed C.bed -z -c hg38 -V 3 
            --force-chrom-peak --force-chrom-more-bed

  # TODO USE SIMPLE.07 instead !!

This command line will compute intersections of all files in more-bed with the file in input (-p) as if the more-bed were regions specified in a GTF

In this case, it will compute the pairwise enrichment of query with A, wuery with B, and query with C.

I MUST PUT A MORE-BED EXAMPLE !!!!! with simple_07 maybe
RQ : NOW I NO LONGER NEED TO SPECIFY MORE-BED-LABELS NO ? NEED TO AMEND THE DOCUMENTAION AND THE FUNCTION NOTES
TO REFLECT THAT















ologram (multiple overlaps)
~~~~~~~~~~~~~~~~~~~~~~






It is also possible to use the **OLOGRAM-MODL** Multiple Overlap Dictionary Learning) plugin to find multiple overlaps (ie. between n>= 2 sets) enrichment.
This is done on the BEDs supplied with the `--more-bed` argument. 






You can ask for all combinations, but 2**N can be big. 



We also give the option to use sparse dictionary learning on the true overlaps
to identify interesting combinations, but you can also specify them yourself.
Mon algorithme MODL (Multiple Overlap Dictionary Learning) de détection des combinaisons via factorisation matricielle et filtrage par algorithme glouton y est intégré. Pour rappel, cet algorithme ne sert qu'à filtrer l'output d'OLOGRAM en termes de combinaisons affichées (OLOGRAM ne calculera l'enrichissement que des combinaisons jugées intéressantes). Par défaut le programme ne l'utilise pas et renvoie toutes les combinaisons... RENCONTREES DANS LES VRAIES DATA, pas dans les shuffles.
  The parameter to use it is --multiple-overlap-max-combinations

Ceci dit je pense que vous n'en aurez pas trop besoin de MODL dans la plupart des cas.
This is mostly useful if there are many files to reduce the number of displayed combinations.
Unlike classical association rules mining algorithms, this focuses on mining complexes and correlation groups (item sets).
Donc (maybe in ologram.py __notes__ only) Quand vous demander à MODL de restreindre le nombre de combinaisons, demandez le top 20 ou 30 pas plus. C'est fait pour trouver des complexes, pas des règles d'association : si vous demandez plus de combi le temps de calcul augmente de manière exponentielle (heures ou jours !). Si vous les voulez toutes, ignorez MODL.
The idea is to use this algorithm to not have all 2**N combinations show. It is designed to find relevant bio clusters.


Add direct link to dict_learning.py source on github

ADD DETAILS !

Say this :'
  Details are available in the code and paper. Broadly speaking, this algorithm will perform many matrix factorizations on the 
  matrix of true overlaps to identify relevant groups of TRs.
  Then a greedy algorithm based on how much these words improve the reconstruction will select the utmost best words
'


SAY you should not ask MODL for more than 20-50 combinations, it is inefficient with more and not designed for it (improvements pending)

Once interesting combis have been found, we will compute enrichment using the OLOGRAM method for the combinations as usual.



You can ask for all combinations. If you want, we have also added a plugin to not show all 2^N combinations (for N files).
It is done with DL (or apriori now as an option? NO DO NOT USE APRIORI !!!!!) iT IS OPTINAL
  Must say that it is about passing -1 (default) to a parameter max_multi_overlap_combis or something like that
Each combination is of the form A+B+C where A, B and C are bed files given as more-bed. They will each have a p value and NB enrichment.

Acknowledge that this plugin of itemset mining is WIP, but it is only used to display only certain combis (use the word "display")
NEW : you can also use apriori for this purpose with the argument --use-apriori-or-something



I heartily recommend using --bed-incl or --bed-excl to restrict the shuffles (ie. shuffling on enhancers only), otherwise longer combis are statitically very improbable




To use MODL, use the --multiple-overlap-max-number-of-combinations argument, with the wanted number of combinations
Also explain rile of --multiple_overlap_target_combi_size : combis longer than this will be ignored. Useful for exact.





**Exact combinations **: Here explain exact and the three cases (see Zim)
  Actually two, simple ! By default you have inexact combis, meaning that at a given position overlaps of A+B+C will count as one towards A+B+...
  To get eact overlaps (A+B but NOT C), put --target--combi size equal to number of --more-beds plus 1 for the query (in the example above, it would be XXX)
  You will know combis are inexact when ther are "..." in the labels.

In most cases use the -z or --no-gtf argument and only pass --more-bed

**Example:**

.. command-output:: gtftk ologram -z -p simple_07_peaks.bed -c simple_07.chromInfo -u 2 -d 2 -K ologram_output --no-date -k 8 --more-bed simple_07_peaks.1.bed simple_07_peaks.2.bed --more-bed-labels One,Two --more-bed-multiple-overlap
	:shell:



MINIBATCH_NB=10
MINIBATCH_SIZE=100
THREADS=8
QUERY=./source.bed
DATA_FILES_DIR=./data
# Query is the file to compare against. Intersections not including the query file will be discarded
# Data files dir is the path to the directory containing the regions of interest as bed files (A.bed, B.bed, C.bed, etc.)
# The program will return the enrichment of relevant combinations such as Query+A, Query+B+C, etc.
# Run OLOGRAM-MODL
gtftk ologram -z -c hg38 -p ${QUERY} --more-bed `ls -d ${DATA_FILES_DIR}/*` 
  -o results --force-chrom-peak --force-chrom-more-bed 
  -V 3 -k ${THREADS} -mn ${MINIBATCH_NB} -ms ${MINIBATCH_SIZE} --more-bed-multiple-overlap
# To use the MODL combination filtering algorithm, add the --multiple-overlap-max-number-of-combinations 42 argument to the previous command line, replacing 42 with the wanted number of combinations
# Also explain rile of --multiple_overlap_target_combi_size

ADD AN EXAMPLE WITH the --multiple-overlap-max-combinations-or-something ARGUMENT
AND ANOTHER EXAMPKE WITH THE --max-size-of-combi-or-something ARGUMENT TO EXPLAIN EXACT AND INEXACT !!!
  gtftk ologram -z -c hg38 -p {input} |\                      # The query
      --more-bed {params.trs} 
      -o results --force-chrom-peak --force-chrom-more-bed  |\
      -V 3 -k 8 -mn 40 -ms 10 |\          # Verbosity, threads, number and size of minibatches
      --more-bed-multiple-overlap         # Take multiple overlaps
      --multiple-overlap-max-number-of-combinations 10     # OPTIONAL ARGUMENT. Use MODL to restrict to THIS MANY combinations (optional)
      --multiple-overlap-target-combi-size 3               # OPTIONAL ARGUMENT. Combis restricted to this size. Also Explain exact (optional)



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


  NOTE : I ONLY NEED TO SHOW ONE QUICK EXAMPLE.
    JUST SHOW THE OLOGRAM RESULT HERE on simple_07










MODL can also be used independantly as a combination mining algorithm.

You need data with one line per transaction and one column per element

For more details, see code comments and paper.


SHOULD THIS GO IN API.RST INSTEAD ?

``` python
    >>> from pygtftk.stats.intersect.dict_learning import Modl, test_data_for_modl
    >>> import numpy as np
    >>> np.random.seed(42)
    >>> flags_matrix = test_data_for_modl(nflags = 1000, number_of_sets = 6, noise = 0.1, cor_groups = [(0,1),(0,1,2,3),(4,5)])
    >>> combi_miner = Modl(flags_matrix, 
    >>>        multiple_overlap_target_combi_size = -1,    # Limit the size of the combinations
    >>>        multiple_overlap_max_number_of_combinations = 3,    # How many words to find ?
    >>>        nb_threads = 1,
    >>>        step_1_factor_allowance = 2)    # How many words to ask for in each step 1 rebuilding
    >>> interesting_combis = combi_miner.find_interesting_combinations()
    >>> assert set(interesting_combis) == set([(1,1,0,0,0,0),(1,1,1,1,0,0),(0,0,0,0,1,1)])
    
```





Please read the nodes below for more details !




**Arguments:**

.. command-output:: gtftk ologram -h
	:shell:













WARNING : if using lots of file, modl may clog and have too big of a matrix !!
Then you should specify custom combis only (show how)

































ologram_merge_stats
~~~~~~~~~~~~~~~~~~~~~~

**Description:** Merge results from different *OLOGRAM* calls in a heatmap for visualisation.


Can still work with OLOGRAM-MODL type results, since they follow the same basic format of one element/combination per line.



.. command-output:: gtftk ologram_merge_stats H3K4me3_ologram_stats.tsv H3K36me3_ologram_stats.tsv H3K79me2_ologram_stats.tsv -o merged_ologram.pdf --labels H3K4me3,H3K36me3,H3K79me2
	:shell:


.. raw:: html

  <br>
  <table>
  <tr>
  <td valign="top">
  <iframe src="_static/example_pa_05.pdf" title="your_title" align="top" width="500" height="620" width="50%" frameborder="0" scrolling="auto" target="Message">
  </iframe>
  </td>
  </tr>
  </table>
  <br>
  <br>

This also works on multiple overlap results

**Arguments:**

.. command-output:: gtftk ologram_merge_stats -h
	:shell:











ologram_modl_treeify
~~~~~~~~~~~~~~~~~~~~~~

**Description:** visualize n-wise enrichment results as a tree by showing strength of association between sets (based on S p-val). Sort of a correlation network.
Will also give a tree of combinations.

Hmm now it gives only the tree of combinations, that first tree is actually garbage I think.


Works on the result (tsv file) of an ologram call with --multiple-overlap


SHOW THE RESULT HERE QUICKLY ON SIMPLE_07

label is optional


SAY IT IS THE PREFERRED REPRESENTATINO FOR OLOGRAM multiple overlap results


Remember that you can EDIT the tsv before passing it to ologram_modl_treeify, for example keeping only the combinations you want


.. command-output:: gtftk ologram_merge_stats -h
	:shell:
# Grab newest tsv file and turn it into a tree to visualize the results
gtftk ologram_modl_treeify -i ologram_result.tsv -o ./results/treeified.pdf -l ThisWasTheNameOfTheQuery





SHOW A QUICK EXAMPLE !!!!!

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

Explain : S is total nb of overlapping base pair in reality, fold change is when comapred to shuffle, p value is such


ologram_merge_runs
~~~~~~~~~~~~~~~~~~~~~~

**Description:** to save memory, merge several runs of OLOGRAM into one run, treating each separate run as a super batch of shuffles




OLOGRAM remembers all intersections occuring inside all minibatches to calculate statistics. If you are using
a large number of shuffles and/or very large files, this may cost a lot of RAM.

In practice, you should not need to use more than 500 shuffles. But if you absolutely require increased precision, 
you can run OLOGRAM several times, treat each run as a "batch of batches" and merge and recalculate stats on the merged superbatch
automatically using this command



```bash
# Make several OLOGRAM runs
N_RUNS = 100
for i in {1..$N_RUNS}
do
   ologram
done
# Possible because each run has a different time and will not overwrite the previous results

# Merge those runs
# use ls to get all files in the directory
gtftk ologram_merge_runs --inputfiles `ls ./results/*.tsv` -o ./merged_batches_result.tsv -V 3

# Treeify and other ologram commands can now be called on the resulting tsv

```
Commands from section 'annotation'
------------------------------------

In the example of this section we will need the following example files:

.. command-output:: gtftk get_example -q -d simple -f '*'
	:shell:

.. command-output:: gtftk get_example -q -d mini_real -f '*'
	:shell:

.. command-output:: gtftk get_example -q -d hg38_chr1 -f '*'
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

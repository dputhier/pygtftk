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


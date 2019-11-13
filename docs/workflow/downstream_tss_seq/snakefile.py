import os

RELEASE = ["94"]
SPECIES = ["homo_sapiens"]
SPECIES_SHORT = "hsapiens"

GV="hg38"

WDIR = os.getcwd()
workdir: WDIR

BIO= ["TEC",
      "antisense",
      "lincRNA",
      "miRNA",
      "misc_RNA",
      "non_stop_decay",
      "nonsense_mediated_decay",
      "processed_pseudogene",
      "processed_transcript",
      "protein_coding",
      "retained_intron",
      "sense_intronic",
      "sense_overlapping",
      "snRNA",
      "snoRNA",
      "transcribed_processed_pseudogene",
      "transcribed_unprocessed_pseudogene",
      "unitary_pseudogene",
      "unprocessed_pseudogene"]

rule final:
    input: expand('output/08_weblogo/{gv}/{s}/{r}/{bio}.png', s=SPECIES, r=RELEASE, gv=GV, bio=BIO)

rule get_fasta:
    output: '{gv}.fa'
    shell: '''
    wget http://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.gv}/bigZips/{wildcards.gv}.fa.gz
    gunzip {wildcards.gv}.fa.gz
    '''

rule retrieveGTF_retrieve:
    output: 'output/01_get_gtf/{s}/{r}.gtf'
    shell: """
    gtftk retrieve -V 1 -C -s {wildcards.s} -r {wildcards.r} -cd > {output}
    """

    
rule selectChrom_select_by_regexp:
    input: 'output/01_get_gtf/{s}/{r}.gtf'
    output: 'output/02_select_chrom/{s}/{r}.gtf'
    shell: """
    gtftk select_by_regexp -V 1  -k chrom -r "chr[0-9XY]+" -i {input} -o {output}
    """

rule del_chrMT:
    input: 'output/02_select_chrom/{s}/{r}.gtf'
    output: 'output/03_del_chrMT/{s}/{r}.gtf'
    shell: """
    gtftk select_by_key -i {input} -k seqid -v chrMT -n -o {output}
    """

rule get_biotype:
    input: 'output/03_del_chrMT/{s}/{r}.gtf'
    output: 'output/04_get_biotype/{bio}/{s}/{r}.gtf'
    shell: """
    gtftk select_by_key -i {input} -k transcript_biotype -v {wildcards.bio} -o {output}
    """

    
rule get_tts:
    input: 'output/04_get_biotype/{bio}/{s}/{r}.gtf'
    output: 'output/05_get_tss/{bio}/{s}/{r}.bed'
    shell: """
    gtftk get_5p_3p_coords -i {input} -v  -n gene_name | sort | uniq > {output}
    """

rule getChromosomeSize:
    output: 'output/05b_chrom_size/chrom_size_{s}_{r}.txt'
    params: v=GV
    shell: """
    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
    "select chrom, size from {params.v}.chromInfo" > {output}
    """

rule slop_tts:
    input: bed='output/05_get_tss/{bio}/{s}/{r}.bed', chr='output/05b_chrom_size/chrom_size_{s}_{r}.txt'
    output: 'output/06_slop_tss/{bio}/{s}/{r}.bed'
    shell: """
    slopBed -i {input.bed} -l 50 -r 50 -s -g {input.chr} > {output}
    """

rule fasta_from_bed:
    input: bed='output/06_slop_tss/{bio}/{s}/{r}.bed', fa='{gv}.fa'
    output: 'output/07_fasta_from_bed/{bio}/{s}/{r}_{gv}.fasta'
    shell: """
    fastaFromBed -s  -fi {input.fa} -bed {input.bed} -fo {output} 
    """

rule weblogo:
    input: 'output/07_fasta_from_bed/{bio}/{s}/{r}_{gv}.fasta'
    output: 'output/08_weblogo/{gv}/{s}/{r}/{bio}.png'
    shell: """
    weblogo -s medium -t {wildcards.bio} -W 20  --format PNG  < {input} > {output} 

    """



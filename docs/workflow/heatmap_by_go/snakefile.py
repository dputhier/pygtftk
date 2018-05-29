import os

RELEASE = ["91"]
SPECIES = ["mus_musculus"]
SPECIES_SHORT = "mmusculus"
BIGWIG="ENCFF670CNM"
GV="mm10"

# DNA binding, apoptosis, Actin cytoskeleton,
# Cell surface
GO_ID = ["0003700", "0097194", "0015629", "0009986"]

WDIR = os.getcwd()
workdir: WDIR


rule final:
    input: expand('output/09_profile/profile_{s}_{r}.pdf', s=SPECIES, r=RELEASE, g=GO_ID)

rule retrieveBigwig:
    output: 'output/00_bigwig/bigwig.bw'
    params: bw=BIGWIG
    shell: """
    cd output/00_bigwig
    wget https://www.encodeproject.org/files/{params.bw}/@@download/{params.bw}.bigWig -O bigwig.bw 
    """


rule retrieveGTF_retrieve:
    output: 'output/01_get_gtf/{s}_{r}.gtf'
    shell: """
    gtftk retrieve -V 2 -C -s {wildcards.s} -r {wildcards.r} -cd > {output}
    """

rule selectChrom_select_by_regexp:
    input: 'output/01_get_gtf/{s}_{r}.gtf'
    output: 'output/02_select_chrom/{s}_{r}.gtf'
    shell: """
    gtftk select_by_regexp -V 2  -k chrom -r "chr[0-9XY]+" -i {input} -o {output}
    """

rule oneTranscriptPerGene_random_tx:
    input: 'output/02_select_chrom/{s}_{r}.gtf'
    output: 'output/03_select_random_transcript/{s}_{r}.gtf'
    shell: """
    gtftk random_tx -V 2 -m 1  -i {input} -o {output}
    """

rule selectByGeneOntology_select_by_go:
    input: 'output/03_select_random_transcript/{s}_{r}.gtf'
    output: 'output/04_select_by_go_{g}/{s}_{r}.gtf'
    params: s=SPECIES_SHORT
    shell: """
    gtftk select_by_go -s {params.s} -g {wildcards.g} -i {input} -o {output}
    """

rule retrieveGeneClasses_tabulate:
    input: 'output/04_select_by_go_{g}/{s}_{r}.gtf'
    output: 'output/05_prepare_classes_{g}/{s}_{r}.txt'
    shell: """
    gtftk tabulate -Hun -k transcript_id -i  {input} | awk 'BEGIN{{OFS="\t"}}{{print $1,"{wildcards.g}"}}' > {output}
    """

rule mergeClasses:
    input: expand('output/05_prepare_classes_{g}/{s}_{r}.txt', s=SPECIES, r=RELEASE, g=GO_ID)
    output: 'output/06_merge/merge_class.txt'
    shell: """
    cat {input} > {output}
    """

rule getChromosomeSize:
    output: 'output/07_chrom_size/chrom_size_{s}_{r}.txt'
    params: v=GV
    shell: """
    mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
    "select chrom, size from {params.v}.chromInfo" > {output}
    """



rule computePromoterCoverage_mk_matrix:
    input:  g='output/03_select_random_transcript/{s}_{r}.gtf', \
            bw='output/00_bigwig/bigwig.bw', \
            s='output/07_chrom_size/chrom_size_{s}_{r}.txt'
    output: 'output/08_prepare_profile/{s}_{r}_matrix.zip'
    shell: """
    gtftk mk_matrix -V 2 -i {input.g} -o output/08_prepare_profile/{wildcards.s}_{wildcards.r}_matrix -u 2500 -d 2500 -k 4 -c  {input.s} {input.bw} 
    """


rule computeMetaplot_profile:
    input: z='output/08_prepare_profile/{s}_{r}_matrix.zip', c='output/06_merge/merge_class.txt'
    output: 'output/09_profile/profile_{s}_{r}.pdf'
    shell: """
    gtftk profile  -V 2 -t {input.c} -i {input.z} -o output/09_profile -g tx_classes -f bwig  -if output/09_profile/profile_{wildcards.s}_{wildcards.r}.pdf -pw 8
    """

# c='output/06_merge/merge_class.txt', \
#

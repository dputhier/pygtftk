import os

RELEASE = ["91"]
SPECIES = ["mus_musculus"]
SPECIES_SHORT = "mmusculus"
BIGWIG = ["ENCFF670CNM", "ENCFF629ABG" ,"ENCFF641ALU", "ENCFF716SHQ"]
LABELS = ["H3K4me3", "H3K27ac", "H3K36me3", "H3K4me1"]
LAB_STR= ",".join(LABELS)

GV="hg38"

GENE_BIOTYPE = ["antisense_RNA",
                "lincRNA",
                "protein_coding"]

WDIR = os.getcwd()
workdir: WDIR


rule final:
    input: expand('output/09_profile/profile_{s}_{r}.pdf', s=SPECIES, r=RELEASE, g=GENE_BIOTYPE)

def get_bw(wildcards):
    return wildcards.bw

rule get_bigwig:
    params: bw=get_bw
    output: 'output/00_bigwig/{bw}.bw'
    shell: """
    cd output/00_bigwig
    wget https://www.encodeproject.org/files/{params.bw}/@@download/{params}.bigWig -O {params}.bw 
    """


rule retrieve_gtf:
    output: 'output/01_get_gtf/{s}_{r}.gtf'
    shell: """
    gtftk retrieve -V 1 -C -s {wildcards.s} -r {wildcards.r} -cd > {output}
    """

rule select_chrom_by_regexp:
    input: 'output/01_get_gtf/{s}_{r}.gtf'
    output: 'output/02_select_chrom/{s}_{r}.gtf'
    shell: """
    gtftk select_by_regexp -V 1  -k chrom -r "chr[0-9XY]+" -i {input} -o {output}
    """

rule select_random_tx:
    input: 'output/02_select_chrom/{s}_{r}.gtf'
    output: 'output/03_select_random_transcript/{s}_{r}.gtf'
    shell: """
    gtftk random_tx -V 1 -m 1  -i {input} -o {output}
    """

rule select_by_gene_biotype:
    input: 'output/03_select_random_transcript/{s}_{r}.gtf'
    output: 'output/04_select_by_genebiotype_{g}/{s}_{r}.gtf'
    shell: """
    gtftk select_by_key -v  {wildcards.g} -k gene_biotype -i {input} -o {output}
    """

rule retrieveGeneClasses_tabulate:
    input: 'output/04_select_by_genebiotype_{g}/{s}_{r}.gtf'
    output: 'output/05_prepare_classes_{g}/{s}_{r}.txt'
    shell: """
    gtftk tabulate -Hun -k transcript_id -i  {input} | awk 'BEGIN{{OFS="\t"}}{{print $1,"{wildcards.g}"}}' > {output}
    """

rule mergeClasses:
    input: expand('output/05_prepare_classes_{g}/{s}_{r}.txt', s=SPECIES, r=RELEASE, g=GENE_BIOTYPE)
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
            bw=expand('output/00_bigwig/{bw}.bw', bw=BIGWIG), \
            s='output/07_chrom_size/chrom_size_{s}_{r}.txt'
    output: 'output/08_prepare_profile/{s}_{r}_matrix.zip'
    params: lab=LAB_STR
    shell: """
    gtftk mk_matrix -l {params.lab} -t transcript -V 1 -i {input.g} -o output/08_prepare_profile/{wildcards.s}_{wildcards.r}_matrix -u 2500 -d 2500 -k 4 -c  {input.s} -y {input.bw} 
    """


rule computeMetaplot_profile:
    input: z='output/08_prepare_profile/{s}_{r}_matrix.zip', c='output/06_merge/merge_class.txt'
    output: 'output/09_profile/profile_{s}_{r}.pdf'
    shell: """
    gtftk profile -nm ranging  -V 1 -t {input.c} -i {input.z} -o output/09_profile -g tx_classes -f bwig  -if output/09_profile/profile_{wildcards.s}_{wildcards.r}.pdf -pw 8
    """

# c='output/06_merge/merge_class.txt', \
#

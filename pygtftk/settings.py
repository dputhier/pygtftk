# ---------------------------------------------------------------
# The script for bash completion
# ---------------------------------------------------------------


COMPLETION_SCRIPT = """
# To activate completion of sub-commands
# Add this script to your ~/.bashrc
 _gtftk()
 {{
     local cur;
     COMPREPLY=();
     cur=${{COMP_WORDS[COMP_CWORD]}};
     if [ ${{COMP_CWORD}} == 1 ]; then
     COMPREPLY=( $( compgen -W '{p}' -- ${{cur}}) )
     else
     COMPREPLY=( $( compgen -f -o plusdirs -- ${{cur}} ) )
     fi
 }}
 complete -o nospace -F _gtftk pygtftk
"""


def get_completion_script():
    import pygtftk.cmd_manager
    prg_list = pygtftk.cmd_manager.CmdManager.cmd_obj_list.keys()
    prg_str = " ".join(sorted(prg_list))
    return COMPLETION_SCRIPT.format(p=prg_str)


"""
#---------------------------------------------------------------
# The list of accepted features
#---------------------------------------------------------------

GTF_ACCEPTED_FEATURE_LIST = ["gene", "transcript", "exon", "CDS",
                             "five_prime_utr", "three_prime_utr",
                             "start_codon", "stop_codon"]

GTF_ACCEPTED_FEATURE_CS = ['"{}"'.format(x) for x in GTF_ACCEPTED_FEATURE_LIST]
GTF_ACCEPTED_FEATURE_CS = ",".join(GTF_ACCEPTED_FEATURE_CS)


#---------------------------------------------------------------
# File Format definition (~/.gtftk.ffd)
#---------------------------------------------------------------

FILE_FORMAT_DEF = ## gtftk file format definition file
## >format_name	file_extension	number_of_columns
## COLUMN	name	type	values	default_value	typed_default_value	mandatory
## KEY	key	source	destinations

>GTF	gtf	9
COLUMN	seqid	S	*	"."	"."	true
COLUMN	source	S	*	"."	"."	true
COLUMN	feature	S	{a}	"."	"."	true
COLUMN	start	I	*	"."	-1	true
COLUMN	end	I	*	"."	-1	true
COLUMN	score	F	*	"."	-1.0	false
COLUMN	strand	C	"+","-"	"."	"."	false
COLUMN	phase	I	0,1,2	"."	-1	false
COLUMN	attributes	A	*	"."	"."	false
KEY	attributes("gene_id")	feature("gene")	feature("transcript"),feature("CDS"),feature("exon"),feature("start_codon"),feature("stop_codon"),feature("five_prime_utr"),feature("three_prime_utr")
KEY	attributes("transcript_id")	feature("transcript")	feature("CDS"),feature("exon"),feature("start_codon"),feature("stop_codon"),feature("five_prime_utr"),feature("three_prime_utr")
KEY	attributes("ccds_id")	feature("CDS")	feature("exon"),feature("start_codon"),feature("stop_codon"),feature("five_prime_utr"),feature("three_prime_utr")

>GFF3	gff3	9
COLUMN	seqid	S	*	"."	true
COLUMN	source	S	*	"."	true
COLUMN	type	S	*	"."	true
COLUMN	start	I	*	"."	true
COLUMN	end	I	*	"."	true
.format(a=GTF_ACCEPTED_FEATURE_CS)
"""

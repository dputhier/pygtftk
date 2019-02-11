/*
 * globals.c
 *
 *  Created on: Sep 4, 2018
 *      Author: fafa
 *
 *  The global variables used in the library.
 *  Unfortunately, the mechanism to parse a binary tree (twalk) doesn't allow
 *  the user to pass any parameter to the parsing function (called by twalk).
 *  So, we must use global variables to pass those parameters, as to get the
 *  results of the parsing process.
 */

#include <libgtftk.h>

/*
 * This is a pointer on the column model of a GTF file. It is initialized by
 * the make_columns() function in column.c source file. It is accessed by most
 * of the functions of the library.
 */
COLUMN **column = NULL;

/*
 * A hashtable used to get a column rank from its name.
 */
struct hsearch_data *column_rank;

/*
 * The number of columns, that is set to 9 for GTF files in column.c source
 * file.
 */
int nb_column;

/*
 * a local copy of the GTF_DATA to process. Used mainly in twalk action
 * functions to get GTF data as GTF indexes contain only row numbers.
 */
GTF_DATA *gtf_d = NULL;

/*
 * a table of exon rows used in add_exon_number function. This table is sorted
 * by genomic location and is used to get the rank of each exon.
 */
SORT_ROW *sort_row = NULL;

/*
 * the number of rows in the previous sort_row table
 */
int nb_sort_row;

/*
 * the name of the "number of exons" attribute
 */
char *enf = NULL;

/*
 * the local copies of min and max number of exons values. Used in
 * select_by_number_of_exons function with a twalk.
 */
int min_noe, max_noe;

/*
 * a ROW_LIST to aggregate all the selected rows (their rank) in several
 * commands
 */
ROW_LIST *row_list;

/*
 * the local copies of min and max transcript size values. Used in
 * select_by_transcript_size function with a twalk.
 */
int min_ts, max_ts;

/*
 * index on transcripts used in select_transcript and convert_to_ensembl
 * functions with twalk
 */
INDEX_ID *tid_index = NULL;

/*
 * index on genes used in convert_to_ensembl function with twalk
 */
INDEX_ID *gid_index = NULL;

/*
 * used in convert_to_ensembl function. If set by twalk action function, that
 * means that the added row (transcript or gene) is the first row of the
 * GTF_DATA so the "data[0]" field in the GTF_DATA must be replaced (by gtf_d0)
 */
GTF_ROW *gtf_d0 = NULL;

/*
 * used in get_list function to get results from twalk action function.
 */
TTEXT *vret = NULL;

/*
 * used in select_by_key function to get the rows that contains an attribute
 */
ROW_LIST *all_rows = NULL;

/*
 * used in select_transcript function to look for the specified transcript
 */
ROW_LIST *test_row_list;

/*
 * used in select_transcript function to get the rows related to the
 * searched transcript (most 5p, longest, or shortest)
 */
ROW_LIST **find_row_list;

/*
 * used in convert_to_ensembl function to add "gene" and "transcript" rows in
 * the final GTF_DATA
 */
int nbrow;

/*
 * used in select_transcript function to ask for the shortest, longest or most
 * 5' transcript of genes. Values are defined in libgtftk.h
 */
int tr_type;

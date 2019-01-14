/*
 * libgtftk.h
 *
 *  Created on: Jan 10, 2017
 *      Author: fafa
 *
 *  Header file for the gtftk library.
 *  Contains all the structure definitions and the prototype declarations.
 */

/*
 * If this flag is set, the library output a very huge debug information used
 * to find memory leaks.
 */
//#define GTFTOOLKIT_DEBUG

#ifndef GTFTOOLKIT_GTFTK_SRC_LIB_LIBGTFTK_H_
#define GTFTOOLKIT_GTFTK_SRC_LIB_LIBGTFTK_H_
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <search.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <unistd.h>

/*
 * Debug of memory allocation. Must be linked with libmemory.so shared library.
 * GTFTOOLKIT_DEBUG must be defined.
 */
#ifdef GTFTOOLKIT_DEBUG
extern void *F_calloc(int nb, int size, char *file, const char *func, int line);
extern void *F_realloc(void *ptr, int size, char *file, const char *func, int line);
extern void F_free(void *ptr, char *file, const char *func, int line);
extern void *F_malloc(int size, char *file, const char *func, int line);
#define calloc(nb, size) F_calloc(nb, size, __FILE__, __func__, __LINE__)
#define realloc(ptr, size) F_realloc(ptr, size, __FILE__, __func__, __LINE__)
#define malloc(size) F_malloc(size, __FILE__, __func__, __LINE__)
#define free(p) F_free(p, __FILE__, __func__, __LINE__)
#endif /* GTFTOOLKIT_DEBUG */

/*
 * constants for transcript selection in select_transcript function
 */
#define SHORTEST_TRANSCRIPT 1
#define LONGEST_TRANSCRIPT 2
#define MOST5P_TRANSCRIPT 3

/*
 * type definition for columns and attributes
 */
#define UNSET -100
#define IRREGULAR -2
#define UNAVAILABLE -1	// .
#define UNKNOWN 0 		// ?
#define INTEGER 1
#define FLOAT 2
#define BOOLEAN 3		// 0 or 1
#define CHAR 4
#define STRING 5
#define UNCONSISTENT 100

/*
 * some useful macros used in get_sequences.c
 */
#define MIN(x, y) (x <= y ? x : y)
#define MAX(x, y) (x > y ? x : y)
#define COMPLEMENT(c) (	c == 'A' ? 'T' : \
						c == 'a' ? 't' : \
						c == 'T' ? 'A' : \
						c == 't' ? 'a' : \
						c == 'G' ? 'C' : \
						c == 'g' ? 'c' : \
						c == 'C' ? 'G' : \
						c == 'c' ? 'g' : c)

/*
 * This structure describes the input (i.e. a GTF/BLASTN plain file or gzipped).
 * It is created by the get_gtf_reader function in get_reader.c source file.
 * gzFile or plain file are set depending on the kind of input (gzip or plain).
 * Even if these two elements are exclusive, they are not in an union to be
 * sure to be compatible with the most of native interfaces.
 */
typedef struct TEXTFILE_READER {
	/*
	 * The file name with its path (or "-" for standard input)
	 */
	char *filename;

	/*
	 * A boolean to tell if it is a gzipped file
	 */
	int gz;

	/*
	 * The gzip file descriptor
	 */
	gzFile gzfile;

	/*
	 * The plain file descriptor
	 */
	FILE *plainfile;
} TEXTFILE_READER;

/*
 * The structure that represents an attribute (key/value) in the last column of
 * a GTF file.
 */
typedef struct ATTRIBUTE {
	char *key, *value;
} ATTRIBUTE;

/*
 * A set of ATTRIBUTE
 */
typedef struct ATTRIBUTES {
	ATTRIBUTE *attr;
	int nb;
} ATTRIBUTES;

/*
 * A structure to store a row from a GTF file.
 */
typedef struct GTF_ROW {
	/*
	 * the 8 first fields of a GTF file row
	 */
	char **field;

	/*
	 * the attributes
	 */
	ATTRIBUTES attributes;

	/*
	 * the rank number of the row in the GTF file
	 */
	int rank;

	/*
	 * the link to the next row
	 */
	struct GTF_ROW *next;

} GTF_ROW;

/*
 * This is the structure that holds data in GTF format. It is also the
 * structure used as input/output for most of the functions of the library. To
 * start using the library, one must call the load_GTF() function with a GTF
 * file name in parameter and gets a pointer on a GTF_DATA in return. Then, all
 * the other functions must be called with this GTF_DATA pointer as input.
 * Their result can be another GTF_DATA pointer that can be used as input for
 * another function of the library.
 */
typedef struct GTF_DATA {
	/*
	 * the number of rows
	 */
	int size;

	/*
	 * a table of rows
	 */
	GTF_ROW **data;

	/*
	 * the comments at the beginning of the file, started with "##"
	 */

} GTF_DATA;

/*
 * This structure represents an index on a column, or on an attribute of the
 * last column.
 */
typedef struct INDEX {
	/*
	 * the name of a column (feature, seqid, ...) or of an attribute (gene_id,
	 * transcript_id, ...)
	 */
	char *key;

	/*
	 * the pointer on a binary tree created with tsearch C function in
	 * the index_row function in column.c source file. This tree contains
	 * ROW_LIST elements described later in this file and that contains a
	 * token and the associated list of row numbers (the rows containing the
	 * token as the value of the key (a column name or an attribute name).
	 */
	void *data;

	/*
	 * a reference to the GTF_DATA on which the index has been made
	 */
	GTF_DATA *gtf_data;

	/*
	 * a pointer on the next index
	 */
	struct INDEX *next;

} INDEX;

/*
 * A structure that contains the information about an index in the column model:
 * the column and the rank as an index can contain several indexes. The index
 * can be accessed as: column[index_id.column]->index[index_id.index_rank]
 */
typedef struct INDEX_ID {
	int column;
	int index_rank;
} INDEX_ID;

/*
 * This is a structure that modelize a column of a GTF file.
 */
typedef struct COLUMN {
	/*
	 * the rank number of the column
	 */
	int num;

	/*
	 * the column name : seqid, source, feature, start, end, score, strand,
	 * phase or attributes
	 */
	char *name;

	/*
	 * the default value to print if no value is available (".")
	 */
	char *default_value;

	/*
	 * a linked list of indexes. It contains only one pointer for each column except
	 * the attributes column for which there is as needed indexes (one can
	 * index data on several attributes)
	 */
	INDEX **index;

	/*
	 * the number of indexes in the previous table
	 */
	int nb_index;
} COLUMN ;

typedef struct STRING_TO_INT_HASH {
	char *key;
	int value;
} STRING_TO_INT_HASH;

/*
 * A list of row numbers associated with a token (the values in the 8 first
 * columns or the values associated to an attribute in the last column). This
 * structure is used in the indexes as elements.
 */
typedef struct ROW_LIST {
	/*
	 * the token that is contained in the rows. For example, this can be "gene"
	 * or "transcript" for an index on the column feature, or "protein_coding"
	 * and "lincRNA" for an index on the attribute "gene_biotype".
	 */
	char *token;

	/*
	 * the number of rows
	 */
	int nb_row;

	/*
	 * the table of row numbers
	 */
	int *row;
} ROW_LIST;

/*
 * This is a structure that can hold any tabulated text. It is for example the
 * result of extract_data function. All functions that return a RAW_DATA
 * structure are "terminal" functions because this kind of result cannot be the
 * input of another function.
 */
typedef struct RAW_DATA {
	/*
	 * The number of rows and columns
	 */
	int nb_rows, nb_columns;

	/*
	 * The name of the columns
	 */
	char **column_name;

	/*
	 * The data (nb_rows x nb_columns character strings)
	 */
	char ***data;
} RAW_DATA;

/*
 * This structure is used to store a list of strings. Useful as a function
 * return type like get_attribute_list in get_list.c source file. Also used as
 * hashtable elements to discard redundant rows in extract_data.
 */
typedef struct STRING_LIST {
	/*
	 * the strings
	 */
	char **list;

	/*
	 * the size of the previous list
	 */
	int nb;
} STRING_LIST;

/*
 * Used by get_sequences function to modelize exons, introns ...
 */
typedef struct SEQFRAG {
	int start, end;
	char strand;
} SEQFRAG ;

typedef struct FEATURE {
	char *name;
	int start, end, tr_start, tr_end;
} FEATURE;

typedef struct FEATURES {
	FEATURE **feature;
	int nb;
} FEATURES;

typedef struct SEQUENCE {
	char *header, *sequence, strand, *seqid, *gene_id, *transcript_id, *gene_name, *gene_biotype;
	int start, end;
	FEATURES *features;
} SEQUENCE;

typedef struct SEQUENCES {
	int nb;
	SEQUENCE **sequence;
} SEQUENCES;

/*
 * Used by get_list function to store an return the results as a matrix of
 * character strings.
 */
typedef struct TTEXT {
	int size;
	char ***data;
} TTEXT;

/*
 * used by add_exon_number to sort exons by their start value
 */
typedef struct SORT_ROW {
	int row;
	int value;
} SORT_ROW;

typedef struct BLAST_HEADER {
	char *program_name;
	char *database_name;
	unsigned int database_length;
	int database_nb_sequences;
} BLAST_HEADER;

typedef struct BLAST_QUERY {
	char *query_name;
	int query_length;
	int nb_subject;
} BLAST_QUERY;

typedef struct BLAST_SUBJECT {
	char *subject_name;
	int subject_length;
	int nb_HSP;
} BLAST_SUBJECT;

typedef struct BLAST_HSP {
	BLAST_HEADER bh;
	BLAST_QUERY bq;
	BLAST_SUBJECT bs;
	double score;
	double expect;
	char *identities;
	int identities_percent;
	char *gaps;
	int gap_percent;
	char strand_query, strand_subject;
	int query_start, query_end, subject_start, subject_end;
} BLAST_HSP;

/*
 * Prototypes for the visible functions (callable by external client)
 */
GTF_DATA *load_GTF(char *input);
GTF_DATA *select_by_key(GTF_DATA *gtf_data, char *key, char *value, int not);
void print_gtf_data(GTF_DATA *gtf_data, char *output, int add_chr);
GTF_DATA *select_by_transcript_size(GTF_DATA *gtf_data, int min, int max);
GTF_DATA *select_by_number_of_exon(GTF_DATA *gtf_data, int min, int max);
GTF_DATA *select_by_genomic_location(GTF_DATA *gtf_data, int nb_loc, char **chr, int *begin_gl, int *end_gl);
RAW_DATA *extract_data(GTF_DATA *gtf_data, char *key, int base, int uniq);
void print_raw_data(RAW_DATA *raw_data, char delim, char *output);
GTF_DATA *select_transcript(GTF_DATA *gtf_data, int type);
SEQUENCES *get_sequences(GTF_DATA *gtf_data, char *genome_file, int intron, int rc);
int free_gtf_data(GTF_DATA *gtf_data);
int free_raw_data(RAW_DATA *raw_data);
char *get_memory(long int size);
int free_mem(char *ptr);
TTEXT *get_feature_list(GTF_DATA *gtf_data);
TTEXT *get_seqid_list(GTF_DATA *gtf_data);
TTEXT *get_attribute_list(GTF_DATA *gtf_data);
TTEXT *get_attribute_values_list(GTF_DATA *gtf_data, char *attribute);
int get_type(GTF_DATA *gtf_data, char *key, int ignore_undef);
GTF_DATA *convert_to_ensembl(GTF_DATA *gtf_data);
GTF_DATA *add_attributes(GTF_DATA *gtf_data, char *features, char *key, char *new_key, char *inputfile_name);
GTF_DATA *del_attributes(GTF_DATA *gtf_data, char *features, char *keys);
GTF_DATA *select_by_positions(GTF_DATA *gtf_data, int *pos, int size);
GTF_DATA *add_exon_number(GTF_DATA *gtf_data, char *exon_number_field);
GTF_DATA *add_prefix(GTF_DATA *gtf_data, char *features, char *key, char *txt, int suffix);
GTF_DATA *merge_attr(GTF_DATA *gtf_data, char *features, char *keys, char *dest_key, char *sep);
GTF_DATA *load_blast(char *input);
GTF_DATA *add_attr_to_pos(GTF_DATA *gtf_data, char *inputfile_name, char *new_key);
void clear_indexes(void);
GTF_DATA *add_attr_column(GTF_DATA *gtf_data, char *inputfile_name, char *new_key);
int int_array_test(int *pos, int size);
void *print_bed(GTF_DATA *gtf_data, char *output, int add_chr, char *keys, char *sep, char *more_info);

#endif /* GTFTOOLKIT_GTFTK_SRC_LIB_LIBGTFTK_H_ */

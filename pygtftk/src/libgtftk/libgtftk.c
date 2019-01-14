/*
 * gtftk_lib.c
 *
 *  Created on: Apr 1, 2016
 *      Author: fafa
 *
 *  The aim of this library is to allow clients written in other languages
 *  than C to take advantage of the speed of C language to access and
 *  manipulate data stored in GTF files. Thus, libgtftk is a library that
 *  provides a set of basic functions to parse, index and extract data from a
 *  GTF file. Those functions can be called from a client written in any
 *  language (C, C++, Java, Python ...) provided that a dedicated C interface
 *  is available (ex : JNI for Java or ctypes for Python). Obviously, C and
 *  C++ clients doesn't need such an interface. Each GTF related function is
 *  implemented in a separate source file. So, this file is just here to hold
 *  common utility functions and most of them are not callable from the shared
 *  library.
 */

#include "libgtftk.h"

extern void print_row(FILE *output, GTF_ROW *r, char delim, int add_chr);

/*
 * This function splits a character string (s) into a table of words (*tab),
 * according to a set of delimiters (delim). Each character of delim is a
 * delimiter. The string is splitted in place and the resulting word table
 * just contains pointers to the words.
 * Parameters:
 * 		tab:	the address of a pointer on a table of string characters
 * 				the table must NOT be reserved before calling this function
 * 		s:		the character string to be splitted
 * 		delim:	the set of charater delimiters
 *
 * Return:		the number of words that were found
 */
int split_ip(char ***tab, char *s, char *delim) {
	int i, n, k, in_token, l;

	in_token = n = k = 0;
	l = strlen(s);
	for (i = 0; i < l; i++)
		if (strchr(delim, (int)(*(s + i))) != NULL) {
			*(s + i) = 0;
			in_token = 0;
		}
		else if (!in_token) {
			in_token = 1;
			n++;
		}
	*tab = (char **)calloc(n, sizeof(char *));
	for (i = 0; i < l; i++)
		if (*(s + i ) != 0) {
			(*tab)[k++] = s + i;
			i += strlen(s + i);
		}
	return n;
}

/*
 * This function trim the extra space characters at the beginning and at the
 * end of the string s. The trimming is performed in place, so no extra memory
 * is reserved to store the trimmed string. The returned value is a pointer on
 * the first non space character of s, or NULL if s contains only space
 * characters.
 *
 * Parameters:
 * 		s:		the string to be trimmed
 *
 * Return:		a pointer on the trimmed string
 */
char *trim_ip(char *s) {
	int b, e, l;

	l = strlen(s);
	for (b = 0; b < l; b++)
		if (*(s + b) != ' ' && *(s + b) != '\t')
			break;
	for (e = l - 1; e > 0; e--)
		if (*(s + e) == ' ' || *(s + e) == '\t')
			*(s + e) = 0;
		else
			break;
	return s + b;
}

/*
 * This function split a key/value pair (as specified in the GTF format for
 * the attributes column) into 2 character strings. value can be delimitted by
 * double quote characters like in Ensembl format or not. The strings pointed
 * by key and value are allocated in this function.
 *
 * Parameters:
 * 		s:		the key/value pair from a GTF attribute
 * 		key:	a pointer to store the address of the key
 * 		value:	a pointer to store the address of the value
 */
void split_key_value(char *s, char **key, char **value) {
	int k = 0;

	if (*s != 0) {
		while (*s == ' ') s++;
		while (*(s + k) != ' ') k++;
		*(s + k) = 0;
		*key = strdup(s);
		s += k + 1;
		while ((*s == ' ') || (*s == '"')) s++;
		k = 0;
		while ((*(s + k) != '"') && (*(s + k) != 0)) k++;
		*(s + k) = 0;
		*value = strdup(s);
	}
}

/*
 * This function is used by the C tsearch and tfind functions to search and
 * add ROW_LIST elements in all the indexes on the columns of a GTF file. For
 * more information, see the man pages of tsearch.
 */
int compare_row_list(const void *p1, const void *p2) {
	ROW_LIST *rl1 = ((ROW_LIST *)p1);
	ROW_LIST *rl2 = ((ROW_LIST *)p2);
	if (!strcmp(rl1->token, "*") || !strcmp(rl2->token, "*")) return 0;
	return strcmp(rl1->token, rl2->token);
}

/*
 * This function is used by the C tsearch function to search and add STRING_LIST
 * elements in a hashtable for removing duplicates rows in extract_data., char *name
 */
int compare_string_list(const void *p1, const void *p2) {
	STRING_LIST *sl1 = ((STRING_LIST *)p1);
	STRING_LIST *sl2 = ((STRING_LIST *)p2);
	int i;

	if (sl1->nb == sl2->nb) {
		for (i = 0; i < sl1->nb; i++)
			if (strcmp(sl1->list[i], sl2->list[i]))
				return strcmp(sl1->list[i], sl2->list[i]);
		return 0;
	}
	return sl1->nb - sl2->nb;
}

/*
 * This function is used by the C qsort function to sort a table of integers
 * representing rows in a GTF file. As some operations can modify the order of
 * the rows in the results, this function is used to be sure to output the
 * rows in the original order.
 */
int comprow(const void *m1, const void *m2) {
	 int *r1 = (int *)m1;
	 int *r2 = (int *)m2;
	 return *r1 - *r2;
}

/*
 * These 2 functions merge two lists of rows (ROW_LIST structures) with
 * suppression of duplications. The resulting list is not sorted.
 *
 * Parameters:
 * 		scr: 	pointer on the source list
 * 		dsr:	pointer on the destination list
 *
 * Return:		the number of rows in dest list, after merging
 */
int add_row(int src, ROW_LIST *dst) {
	if (dst == NULL) {
		dst = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));
		dst->row = (int *)calloc(1, sizeof(int));
		dst->row[dst->nb_row] = src;
		dst->nb_row++;
	}
	else if (bsearch(&src, dst->row, dst->nb_row, sizeof(int), comprow) == NULL) {
		dst->row = (int *)realloc(dst->row, (dst->nb_row + 1) * sizeof(int));
		dst->row[dst->nb_row] = src;
		dst->nb_row++;
	}
	return dst->nb_row;
}

int add_row_list(ROW_LIST *src, ROW_LIST *dst) {
	int i;
	for (i = 0; i < src->nb_row; i++) add_row(src->row[i], dst);
	return dst->nb_row;
}

/*
 * This function prints the content of a GTF_DATA on the standard output. It
 * can be called from the library (default visibility).
 *
 * Parameters:
 * 		gtf_data:	a pointer on the GTF data to be printed
 */
__attribute__ ((visibility ("default")))
void print_gtf_data(GTF_DATA *gtf_data, char *output, int add_chr) {
	int i;
	FILE *out = stdout;

	if (gtf_data != NULL) {
		if (*output != '-') out = fopen(output, "w");
		if (out == NULL) out = stdout;
		for (i = 0; i < gtf_data->size; i++) print_row(out, gtf_data->data[i], '\t', add_chr);
		if (out != stdout) {
			fflush(out);
			fclose(out);
		}
	}
}

/*
 * This function prints the content of a RAW_DATA (a table of string elements)
 * with the given delimiter.
 *
 * Parameters:
 * 		raw_data:	a pointer on the RAW_DATA structure to be printed
 * 		delim:		a single character delimiter
 * 		output:		a file name or "-" for standard output
 */
__attribute__ ((visibility ("default")))
void print_raw_data(RAW_DATA *raw_data, char delim, char *output) {
	int i, k;
	FILE *out = stdout;

	if (*output != '-') out = fopen(output, "w");
	if (out == NULL) out = stdout;
	fprintf(out, "%s", raw_data->column_name[0]);
	for (i = 1; i < raw_data->nb_columns; i++) fprintf(out, "%c%s", delim, raw_data->column_name[i]);
	fprintf(out, "\n");
	for (i = 0; i < raw_data->nb_rows; i++) {
		fprintf(out, "%s", raw_data->data[i][0]);
		for (k = 1; k < raw_data->nb_columns; k++) fprintf(out, "%c%s", delim, raw_data->data[i][k]);
			fprintf(out, "\n");
	}
	if (out != stdout) {
		fflush(out);
		fclose(out);
	}
}

/*
 * Look for an attribute in a row.
 *
 * Parameters:
 * 		row:	the row to look in
 * 		at:		the attribute name to look for
 *
 * Returns:		the rank of the found attribute (or -1 if not found)
 */
int is_in_attrs(GTF_ROW *row, char *at) {
	int ret = -1, i;
	for (i = 0; i < row->attributes.nb; i++)
		if (!strcmp((row->attributes.attr + i)->key, at)) {
			ret = i;
			break;
		}
	return ret;
}

/*
 * Get the value of an attribute in a row
 *
 * Parameters:
 * 		row:	the row in wich to look for the attribute
 * 		attr:	the name of the attribute
 *
 * Returns:		the value of the attribute or NULL if it does not exist
 */
char *get_attribute_value(GTF_ROW *row, char *attr) {
	int k = is_in_attrs(row, attr);

	if (k != -1) return (row->attributes.attr + k)->value;
	return NULL;
}

/*
 * Some convenient methods for a client (in Python for example) to allocate
 * and free memory blocks
 */
__attribute__ ((visibility ("default")))
char *get_memory(long int size) {
	int i;
	char *mem = calloc(size, 1);
	for (i = 0; i < size; i++) mem[i] = i & 0xFF;
	return mem;
}

__attribute__ ((visibility ("default")))
int free_mem(char *ptr) {
	free(ptr);
	return 0;
}

/*
 * This function make a clone of the input GTF data. Used in functions that
 * does annotation tasks (add_attributes, add_exon_number ...). Those functions
 * make a clone of the input GTF data and perform the annotation on the clone,
 * leaving the initial data unchanged.
 */
GTF_DATA *clone_gtf_data(GTF_DATA *gtf_data) {
	int i, j;
	GTF_ROW *row;

	/*
	 * The clone GTF data to return
	 */
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));
	ret->size = gtf_data->size;
	ret->data = (GTF_ROW **)calloc(ret->size, sizeof(GTF_ROW *));

	/*
	 * The GTF rows loop
	 */
	for (i = 0; i < ret->size; i++) {
		/*
		 * creates a new row, put it in the cloned gtf_data, and fill it
		 */
		row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));
		ret->data[i] = row;
		row->rank = gtf_data->data[i]->rank;
		row->attributes.nb = gtf_data->data[i]->attributes.nb;

		/*
		 * make the rows linked list at the same time
		 */
		if (i > 0) ret->data[i - 1]->next = row;

		/*
		 * creates the attributes table
		 */
		//row->attributes.attr = (ATTRIBUTE **)calloc(row->attributes.nb, sizeof(ATTRIBUTE *));
		row->attributes.attr = (ATTRIBUTE *)calloc(row->attributes.nb, sizeof(ATTRIBUTE));
		//row->attributes.attr = (ATTRIBUTE **)bookmem(row->attributes.nb, sizeof(ATTRIBUTE *), __FILE__, __func__, __LINE__);

		/*
		 * the attributes loop
		 */
		for (j = 0; j < row->attributes.nb; j++) {
			/*
			 * create an attribute in the attributes table of the row and fill it
			 */
			//row->attributes.attr[j] = (ATTRIBUTE *)calloc(1, sizeof(ATTRIBUTE));
			//row->attributes.attr[j]->value = strdup(gtf_data->data[i]->attributes.attr[j]->value);
			//row->attributes.attr[j]->key = strdup(gtf_data->data[i]->attributes.attr[j]->key);
			(row->attributes.attr + j)->value = strdup((gtf_data->data[i]->attributes.attr + j)->value);
			(row->attributes.attr + j)->key = strdup((gtf_data->data[i]->attributes.attr + j)->key);

			/*
			 * make the attributes linked list at the same time
			 */
			//if (j > 0) row->attributes.attr[j - 1]->next = row->attributes.attr[j];
		}
		/*
		 * copy the 8 first columns values
		 */
		row->field = (char **)calloc(8, sizeof(char *));
		for (j = 0; j < 8; j++) row->field[j] = strdup(gtf_data->data[i]->field[j]);
	}

	return ret;
}
/*
 * Add an attribute in a row. The memory needed to store the attribute is
 * reserved and the new attribute is linked with the last one.
 */
void add_attribute(GTF_ROW *row, char *key, char *value) {
	row->attributes.nb++;
	/*row->attributes.attr = (ATTRIBUTE **)realloc(row->attributes.attr, row->attributes.nb * sizeof(ATTRIBUTE *));
	row->attributes.attr[row->attributes.nb - 1] = (ATTRIBUTE *)calloc(1, sizeof(ATTRIBUTE));
	row->attributes.attr[row->attributes.nb - 1]->key = strdup(key);
	row->attributes.attr[row->attributes.nb - 1]->value = strdup(value);
	if (row->attributes.nb > 1)
		row->attributes.attr[row->attributes.nb - 2]->next = row->attributes.attr[row->attributes.nb - 1];*/
	row->attributes.attr = (ATTRIBUTE *)realloc(row->attributes.attr, row->attributes.nb * sizeof(ATTRIBUTE));
	(row->attributes.attr + row->attributes.nb - 1)->key = strdup(key);
	(row->attributes.attr + row->attributes.nb - 1)->value = strdup(value);
}

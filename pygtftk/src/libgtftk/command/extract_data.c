/*
 * extract_data.c
 *
 *  Created on: Mar 31, 2016
 *      Author: fafa
 *
 * The extract_data function returns a RAW_DATA structure that contains a
 * selection of information from a GTF_DATA structure (some columns or some
 * attributes).
 */

#include  "libgtftk.h"

/*
 * external functions declaration
 */
extern int split_ip(char ***tab, char *s, char *delim);
extern STRING_LIST *get_all_attributes(GTF_DATA *gtf_data);
extern int is_in_attrs(GTF_ROW *row, char *at);
extern int compare_string_list(const void *p1, const void *p2);

/*
 * global variables declaration
 */
extern COLUMN **column;
extern int nb_column;

/*
 * This function is called by twalk on the hashtable used to remove redundant
 * rows in the result when uniq parameter is 1. As the tdestroy function is not
 * available in OS X, we must do the job !!!
 */
static void destroy_string_list_tree(const void *nodep, const VISIT which, const int depth) {
	STRING_LIST *sl = *(STRING_LIST **)nodep;
	int i;

	switch (which) {
		case preorder:
			break;
		case postorder:
			break;
		case endorder:
		case leaf:
			for (i = 0; i < sl->nb; i++) free(sl->list[i]);
			free(sl->list);
			free(sl);
			break;
	}
}

/*
 * Look for a column in the column names.
 *
 * Parameters:
 * 		col:	the column name to look for
 *
 * Returns:		the rank of the found column (or -1 if not found)
 */
int is_in_columns(char *col) {
	int ret = -1, i;
	for (i = 0; i < nb_column; i++)
		if (!strcmp(column[i]->name, col)) {
			ret = i;
			break;
		}
	return ret;
}

/*
 * This function extracts data from a GTF_DATA structure and pack it into
 * a RAW_DATA structure, that is a matrix of strings.
 *
 * Parameters:
 * 		gtf_data:	the GTF source data
 * 		key:		the list of column names and attribute names that are to be
 * 					extracted (separated by "," character)
 * 		base:		0 or 1 : 0-based or 1-based data
 * 		uniq:		0 or 1 : keep or discard redundant rows in results
 *
 * Returns:			a matrix of strings, packed into a RAW_DATA structure
 */
__attribute__ ((visibility ("default")))
RAW_DATA *extract_data(GTF_DATA *gtf_data, char *key, int base, int uniq) {
	/*
	 * Declaration and reservation of the structure to be returned
	 */
	RAW_DATA *ret = (RAW_DATA *)calloc(1, sizeof(RAW_DATA));

	/*
	 * Some convenient local variables
	 */
	int i, k, n, j;

	/*
	 * The list of all attributes in the given GTF_DATA
	 */
	STRING_LIST *attributes;

	/*
	 * hashtable used to discard redundant rows if uniq parameter is set to 1
	 */
	void *hash = NULL;

	/*
	 * a string list to look for i the hashtable
	 */
	STRING_LIST *row, *addedrow;

	/*
	 * if key is "all", we extract ALL the data in the GTF_DATA
	 */
	if (!strcmp(key, "all")) {
		/*
		 * get the list of all attributes in the GTF_DATA
		 */
		attributes = get_all_attributes(gtf_data);

		/*
		 * reserve the memory for the list of column names, in the results
		 * (ROW_DATA)
		 */
		ret->column_name = (char **)malloc((8 + attributes->nb) * sizeof(char *));

		/*
		 * Fill the list of result columns with the column names and the
		 * attribute names of the GTF_DATA
		 */
		for (i = 0; i < 8; i++) ret->column_name[ret->nb_columns++] = column[i]->name;
		for (i = 0; i < attributes->nb; i++) ret->column_name[ret->nb_columns++] = strdup(attributes->list[i]);
		free(attributes->list);
	}
	else {
		/*
		 * if key contains a list of column and attribute names, split this
		 * list right in the RAW_DATA structure
		 */
		ret->nb_columns = split_ip(&(ret->column_name), strdup(key), ",");
	}

	/*
	 * The loop on the rows of GTF_DATA
	 */
	for (k = 0; k < gtf_data->size; k++) {
		/*
		 * reserve the memory for the corresponding result row
		 */
		row = (STRING_LIST *)calloc(1, sizeof(STRING_LIST));
		row->nb = ret->nb_columns;
		row->list = (char **)calloc(row->nb, sizeof(char *));

		/*
		 * The loop on the columns and attributes to extract
		 */
		for (i = 0; i < ret->nb_columns; i++) {

			/*
			 * column extraction
			 */
			if ((n = is_in_columns(ret->column_name[i])) != -1) {
				/*
				 * do a copy of the field, as it may be altered
				 */
				row->list[i] = strdup(gtf_data->data[k]->field[n]);

				/*
				 * if we want data 0-based, we must remove 1 from the start
				 * value
				 */
				if (!strcmp(ret->column_name[i], "start") && !base) {
					j = atoi(row->list[i]) - 1;
					sprintf(row->list[i], "%d", j);
				}
			}

			/*
			 * attribute extraction
			 */
			else if ((n = is_in_attrs(gtf_data->data[k], ret->column_name[i])) != -1)
				//row->list[i] = strdup(gtf_data->data[k]->attributes.attr[n]->value);
				row->list[i] = strdup((gtf_data->data[k]->attributes.attr + n)->value);

			/*
			 * not a column and not an attribute ! (some attributes are not
			 * present in all the rows of a GTF_DATA)
			 */
			else
				row->list[i] = strdup("?");
		}
		/*
		 * if uniq is set, put the row in a hashtable to allow removing of
		 * doubles
		 */
		if (uniq) {
			addedrow = *(STRING_LIST **)tsearch(row, &hash, compare_string_list);
			if (addedrow != row) {
				/*
				 * the row is already present in the results, don't need to add
				 * it again
				 */
				for (i = 0; i < ret->nb_columns; i++) free(row->list[i]);
				free(row->list);
				free(row);
			}
			else {
				/*
				 * reserve the memory for the new row to add and for the
				 * corresponding string list
				 */
				ret->data = (char ***)realloc(ret->data, (ret->nb_rows + 1) * sizeof(char **));
				ret->data[ret->nb_rows] = (char **)calloc(ret->nb_columns, sizeof(char *));
				for (i = 0; i < ret->nb_columns; i++) ret->data[ret->nb_rows][i] = strdup(row->list[i]);
				ret->nb_rows++;
			}
		}
		else {
			/*
			 * reserve the memory for the new row to add and for the
			 * corresponding string list
			 */
			ret->data = (char ***)realloc(ret->data, (ret->nb_rows + 1) * sizeof(char **));
			ret->data[ret->nb_rows] = (char **)calloc(ret->nb_columns, sizeof(char *));
			for (i = 0; i < ret->nb_columns; i++) ret->data[ret->nb_rows][i] = strdup(row->list[i]);
			ret->nb_rows++;
		}
	}

	/*
	 * destroy rows hashtable if needed
	 */
	if (uniq) twalk(hash, destroy_string_list_tree);

	/*
	 * Returns the extracted data
	 */
	return ret;
}

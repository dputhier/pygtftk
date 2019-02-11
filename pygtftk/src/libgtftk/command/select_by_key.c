/*
 * select_by_key.c
 *
 *  Created on: Jan 25, 2017
 *      Author: fafa
 *
 * This is the implementation of select_by_key function.
 * This function selects rows in GTF data that contains (or not) some values
 * for a column or for an attribute.
 */

#include "libgtftk.h"

/*
 * external functions in libgtftk.c
 */
extern int comprow(const void *m1, const void *m2);
extern int compare_row_list(const void *p1, const void *p2);
extern int add_row_list(ROW_LIST *src, ROW_LIST *dst);
extern int split_ip(char ***tab, char *s, char *delim);

/*
 * external functions in load_gtf.c
 */
extern INDEX_ID *index_gtf(GTF_DATA *gtf_data, char *key);
extern int update_row_table(GTF_DATA *gtf_data);

/*
 * global variables
 */
extern COLUMN **column;
extern ROW_LIST *all_rows;

/*
 * This function is used through a twalk in an index on an attribute to get
 * all the rows that really have the attribute. The rows are put in the global
 * variable all_rows.
 */
static void get_all_rows(const void *nodep, const VISIT which, const int depth) {
	ROW_LIST *rl = *(ROW_LIST **)nodep;

	switch (which) {
		case preorder:
			break;
		case postorder:
			break;
		case endorder:
		case leaf:
			add_row_list(rl, all_rows);
			break;
	}
}


/*
 * select_by_key function selects rows in GTF_DATA that contains some given
 * value(s) for a column or an attribute. For instance:
 *  - to get all "gene" and "transcript" rows:
 *  	 select_by_key(gtf_data, "feature", "gene,transcript", 0)
 * 	- to get all the rows concerning the genes ESR1 and ESR2:
 * 		select_by_key(gtf_data, "gene_name", "ESR1,ESR2", 0)
 * 	- to search for all lincRNA genes:
 * 		select_by_key(gtf_data, "gene_biotype", "lincRNA", 0)
 * 	- to get all non coding genes:
 * 		select_by_key(gtf_data, "gene_biotype", "protein_coding", 1)
 * As output from select_by_key can be used as input, one can chain several
 * function calls. For instance, to get the first exon af all transcripts of
 * gene ESR2:
 * 	select_by_key(select_by_key(select_by_key(gtf_data,
 * 		"gene_name", "ESR2", 0),
 * 		"feature", "exon", 0),
 * 		"exon_number", "1", 0)
 *
 * Parameters:
 * 		gtf_data:	a GTF_DATA structure obtained by a call to loadGTF function
 * 		key:		a column name or an attribute name
 * 		value:		a set of values (separated by "," characters) for key to
 * 					look for in the GTF_DATA
 * 		not:		1 to get the rows that don't contains the requested values
 * 					2 to get the rows that don't contains the requested values
 * 					but have the attribute
 *
 * Returns:		a GTF_DATA structure that contains the result of the query
 */
__attribute__ ((visibility ("default")))
GTF_DATA *select_by_key(GTF_DATA *gtf_data, char *key, char *value, int not) {
	int i, j, k, p, n = 0;
	/*
	 * reserve memory for the GTF_DATA structure to return
	 */
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));

	/*
	 * Some ROW_LIST variables needed to perform different tasks:
	 * 	- test_row_list will be used to look for a ROW_LIST element in an index
	 * 	- row_list is used to build the list of rows to keep in result
	 * 	- find_row_list contains the found row list for each value, that will
	 * 	  be merged into row_list
	 */
	ROW_LIST *test_row_list = calloc(1, sizeof(ROW_LIST)), *row_list = NULL, **find_row_list = NULL;

	/*
	 * splitting the given values with the "," character
	 */
	char **values;
	int nb_value = split_ip(&values, value, ",");

	/*
	 * indexing the GTF_DATA with the given key (column name or attribute name)
	 */
	INDEX_ID *index_id = index_gtf(gtf_data, key);

	GTF_ROW *row, *previous_row = NULL;

	/*
	 * allocating memory for the final ROW_LIST
	 */
	row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));

	for (p = 0; p < nb_value; p++) {
		/*
		 * for each value, we look for it in the index returned by index_gtf
		 * and if we find it, we merge (add_row_list) the returned ROW_LIST
		 * (find_row_list) with the final one (row_list). Duplicates are
		 * suppressed by add_row_list function.
		 */
		test_row_list->token = values[p];
		find_row_list = (ROW_LIST **)tfind(test_row_list, &(column[index_id->column]->index[index_id->index_rank]->data), compare_row_list);
		if (find_row_list != NULL) add_row_list(*find_row_list, row_list);
	}

	/*
	 * as searches were made in the order of the given values, row_list may
	 * contains unsorted rows, so we need to sort them in ascending order for
	 * the output not to be mixed
	 */
	qsort(row_list->row, row_list->nb_row, sizeof(int), comprow);

	/*
	 * now it's time to build the GTF_DATA result
	 */
	if (!not) {
		/*
		 * if the query was positive, the size of the result is the number of
		 * rows in row_list
		 */
		ret->size = row_list->nb_row;

		/*
		 * we reserve memory for the first row in the table of rows
		 */
		ret->data = (GTF_ROW **)calloc(1, sizeof(GTF_ROW *));

		/*
		 * the main loop on rows to put in results
		 */
		for (j = 0; j < ret->size; j++) {
			/*
			 * reserve memory for a new row
			 */
			row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));

			/*
			 * at this stage, we store only the first row. The whole table will
			 * be expanded later
			 */
			if (j == 0) ret->data[0] = row;

			/*
			 * set the rank of the row (in the GTF file) and the number of attributes
			 * each row in row_list is a number used to get the real GTF_ROW in the
			 * whole GTF data
			 */
			row->rank = gtf_data->data[row_list->row[j]]->rank;
			row->attributes.nb = gtf_data->data[row_list->row[j]]->attributes.nb;

			/*
			 * reserve memory for the 8 first fields and fill them
			 */
			row->field = (char **)calloc(8, sizeof(char*));
			for (i = 0; i < 8; i++) row->field[i] = strdup(gtf_data->data[row_list->row[j]]->field[i]);

			/*
			 * reserve memory for attributes and fill them
			 */
			/*row->attributes.attr = (ATTRIBUTE **)calloc(row->attributes.nb, sizeof(ATTRIBUTE *));
			for (i = 0; i < row->attributes.nb; i++) {
				row->attributes.attr[i] = (ATTRIBUTE *)calloc(1, sizeof(ATTRIBUTE));
				row->attributes.attr[i]->key = strdup(gtf_data->data[row_list->row[j]]->attributes.attr[i]->key);
				row->attributes.attr[i]->value = strdup(gtf_data->data[row_list->row[j]]->attributes.attr[i]->value);
			}*/
			row->attributes.attr = (ATTRIBUTE *)calloc(row->attributes.nb, sizeof(ATTRIBUTE));
			for (i = 0; i < row->attributes.nb; i++) {
				(row->attributes.attr + i)->key = strdup((gtf_data->data[row_list->row[j]]->attributes.attr + i)->key);
				(row->attributes.attr + i)->value = strdup((gtf_data->data[row_list->row[j]]->attributes.attr + i)->value);
			}

			/*
			 * link the new row with the previous one
			 */
			if (j > 0) previous_row->next = row;
			previous_row = row;
		}
	}
	else {
		/*
		 * the query is negative, so we want to keep all the rows that are NOT
		 * in row_list, then the size of the result is the difference between
		 * the total number of rows and the number of rows in row_list. This
		 * total number of rows depends on the value of the "not" parameter.
		 */
		all_rows = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));
		if (not == 2) {
			/*
			 * We start with all the rows that have the key
			 */
			twalk(column[index_id->column]->index[index_id->index_rank]->data, get_all_rows);
			qsort(all_rows->row, all_rows->nb_row, sizeof(int), comprow);
		}
		else {
			/*
			 * We start with all the rows, even those that don't have the key
			 */
			all_rows->row = (int *)calloc(gtf_data->size, sizeof(int));
			all_rows->nb_row = gtf_data->size;
			for (i = 0; i < all_rows->nb_row; i++) all_rows->row[i] = i;
		}

		/*
		 * compute the size of the results
		 */
		ret->size = all_rows->nb_row - row_list->nb_row;

		/*
		 * we reserve memory for the first row in the table of rows
		 */
		ret->data = (GTF_ROW **)calloc(1, sizeof(GTF_ROW *));

		/*
		 * an ugly code to get the "complement" rows of row_list in gtf_data
		 * should be rewrote asap ...
		 */
		j = 0;
		if (row_list->nb_row == 0) {
			row_list->row = (int *)calloc(1, sizeof(int));
			row_list->row[0] = -1;
		}

		for (k = 0; k < all_rows->nb_row && j < row_list->nb_row; k++) {
			if (all_rows->row[k] < row_list->row[j]) {
				row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));
				if (n == 0) ret->data[0] = row;
				row->field = (char **)calloc(8, sizeof(char*));
				for (i = 0; i < 8; i++) row->field[i] = strdup(gtf_data->data[all_rows->row[k]]->field[i]);
				row->attributes.nb = gtf_data->data[all_rows->row[k]]->attributes.nb;
				/*row->attributes.attr = (ATTRIBUTE **)calloc(row->attributes.nb, sizeof(ATTRIBUTE *));
				for (i = 0; i < row->attributes.nb; i++) {
					row->attributes.attr[i] = (ATTRIBUTE *)calloc(1, sizeof(ATTRIBUTE));
					row->attributes.attr[i]->key = strdup(gtf_data->data[all_rows->row[k]]->attributes.attr[i]->key);
					row->attributes.attr[i]->value = strdup(gtf_data->data[all_rows->row[k]]->attributes.attr[i]->value);
				}*/
				row->attributes.attr = (ATTRIBUTE *)calloc(row->attributes.nb, sizeof(ATTRIBUTE));
				for (i = 0; i < row->attributes.nb; i++) {
					(row->attributes.attr + i)->key = strdup((gtf_data->data[all_rows->row[k]]->attributes.attr + i)->key);
					(row->attributes.attr + i)->value = strdup((gtf_data->data[all_rows->row[k]]->attributes.attr + i)->value);
				}
				row->rank = gtf_data->data[all_rows->row[k]]->rank;
				if (n > 0) previous_row->next = row;
				previous_row = row;
				n++;
			}
			else if (all_rows->row[k] == row_list->row[j])
				j++;
		}
		if (n != ret->size) {
			for (k = row_list->row[row_list->nb_row == 0 ? 0 : row_list->nb_row - 1] + 1; k <= all_rows->row[all_rows->nb_row - 1]; k++) {
				row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));
				if (n == 0) ret->data[0] = row;
				row->field = (char **)calloc(8, sizeof(char*));
				for (i = 0; i < 8; i++) row->field[i] = strdup(gtf_data->data[k]->field[i]);
				row->attributes.nb = gtf_data->data[k]->attributes.nb;
				/*row->attributes.attr = (ATTRIBUTE **)calloc(row->attributes.nb, sizeof(ATTRIBUTE *));
				for (i = 0; i < row->attributes.nb; i++) {
					row->attributes.attr[i] = (ATTRIBUTE *)calloc(1, sizeof(ATTRIBUTE));
					row->attributes.attr[i]->key = strdup(gtf_data->data[k]->attributes.attr[i]->key);
					row->attributes.attr[i]->value = strdup(gtf_data->data[k]->attributes.attr[i]->value);
				}*/
				row->attributes.attr = (ATTRIBUTE *)calloc(row->attributes.nb, sizeof(ATTRIBUTE));
				for (i = 0; i < row->attributes.nb; i++) {
					(row->attributes.attr + i)->key = strdup((gtf_data->data[k]->attributes.attr + i)->key);
					(row->attributes.attr + i)->value = strdup((gtf_data->data[k]->attributes.attr + i)->value);
				}
				row->rank = gtf_data->data[k]->rank;
				if (n > 0) previous_row->next = row;
				previous_row = row;
				n++;
			}
		}
	}

	/*
	 * expanding rows table
	 */
	update_row_table(ret);

	/*
	 * freeing memory
	 */
	free(values);
	free(test_row_list);
	free(row_list->row);
	free(row_list);
	if (all_rows != NULL) {
		if (all_rows->row != NULL) free(all_rows->row);
		free(all_rows);
		all_rows = NULL;
	}
	free(index_id);

	return ret;
}

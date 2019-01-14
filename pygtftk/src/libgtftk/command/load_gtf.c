/*
 * load_gtf.c
 *
 *  Created on: Jan 12, 2017
 *      Author: fafa
 *
 * This source file contains functions to load a GTF file into memory, to get
 * a GTF_DATA structure and to index the whole thing with a given column name
 * or attribute name. The GTF_DATA structure is the object that must be used
 * as input with most of the functions of the library.
 */

#include "libgtftk.h"
#include <search.h>
#include <time.h>

/*
 * external functions in gtf_reader.c
 */
extern TEXTFILE_READER *get_gtf_reader(char *query);
extern char *get_next_gtf_line(TEXTFILE_READER *gr, char *buffer);

/*
 * external functions in column.c
 */
extern void make_columns(void);
extern void make_index(INDEX_ID *index_id, char *key);
extern void index_row(int row_nb, char *value, INDEX *index);
extern int comprow(const void *m1, const void *m2);

/*
 * external functions in libgtftk.c
 */
extern int split_ip(char ***tab, char *s, char *delim);
extern void split_key_value(char *s, char **key, char **value);
extern char *trim_ip(char *);

/*
 * global variables
 */
extern COLUMN **column;
extern int nb_column;

/*
 * Reserve memory for a new attribute index and add it at the end of the
 * attributes column. Used by index_gtf in this source file.
 *
 * Parameters:
 * 		key:	an attribute name
 *
 * Returns:		the rank of the index in the index table of the column
 */
int add_index(char *key, GTF_DATA *gtf_data) {
	column[8]->nb_index++;
	column[8]->index = realloc(column[8]->index, column[8]->nb_index * sizeof(INDEX));
	column[8]->index[column[8]->nb_index - 1]->key = strdup(key);
	column[8]->index[column[8]->nb_index - 1]->data = NULL;
	column[8]->index[column[8]->nb_index - 1]->gtf_data = gtf_data;
	return column[8]->nb_index - 1;
}

/*
 * This comparison function is used in get_all_attributes function to get the
 * list of all the attributes that are in a GTF data. Two attributes are equals
 * if they have the same name.
 */
int compatt(const void *s1, const void *s2) {
	char *cs1 = *(char **)s1;
	char *cs2 = *(char **)s2;
	return strcmp(cs1, cs2);
}

/*
 * In a GTF_DATA, the rows are stored in two manners: a linked list and a table.
 * The linked list is usefull for adding rows (ie convert_to_ensembl function)
 * while the table is quicker to address a particular row.
 * After some operations, the rows can be stored only in the linked list (a
 * pointer on the first row is stored in gtf_data->data[0]). This function
 * rebuild the corresponding table in gtf_data->data from the linked list.
 *
 * Parameters:
 * 		gtf_data:	the input GTF data
 *
 * Returns:			0 if success
 */
int update_row_table(GTF_DATA *gtf_data) {
	GTF_ROW *row = gtf_data->data[0];
	int i;

	gtf_data->data = (GTF_ROW **)realloc(gtf_data->data, gtf_data->size * sizeof(GTF_ROW *));
	for (i = 0; i < gtf_data->size; i++) {
		gtf_data->data[i] = row;
		row = row->next;
	}
	return 0;
}

/*
 * Like the rows, the attributes are stored as a linked list and a table. The
 * linked list is used to remove and add attributes easily, while the table is
 * more suitable to access a single attribute. This function reconstruct the
 * table from the linked list.
 *
 * Parameters:
 * 		row:	a row for which to update the attributes table
 *
 * Returns:		the number of attributes in the row
 */
/*int update_attribute_table(GTF_ROW *row) {
	int i;
	ATTRIBUTE *att;
	att = row->attributes.attr[0];
	free(row->attributes.attr);
	row->attributes.attr = (ATTRIBUTE **)calloc(row->attributes.nb, sizeof(ATTRIBUTE *));
	for (i = 0; i < row->attributes.nb; i++) {
		row->attributes.attr[i] = att;
		att = att->next;
	}
	return row->attributes.nb;
}*/

/*
 * The symetric function of update_row_table. This function rebuild the linked
 * list from the table in gtf_data->data.
 *
 * Parameters:
 * 		gtf_data:	the input GTF data
 *
 * Returns:			0 if success
 */
/*int update_linked_list(GTF_DATA *gtf_data) {
	int i, j;
	GTF_ROW *row;

	for (i = 0; i < gtf_data->size - 1; i++) {
		if ((gtf_data->data[i]->next != NULL) && (gtf_data->data[i]->next != gtf_data->data[i + 1])) {
			row = gtf_data->data[i]->next;
			for (j = 0; j < 8; j++) free(row->field[j]);
			free(row->field);
			for (j = 0; j < row->attributes.nb; j++) {
				free(row->attributes.attr[j]->key);
				free(row->attributes.attr[j]->value);
				free(row->attributes.attr[j]);
			}
			free(row->attributes.attr);
			free(row);
		}
		gtf_data->data[i]->next = gtf_data->data[i + 1];
	}
	return 0;
}*/

/*
 * This function read GTF data from a GTF file (gzipped or not) or standard
 * input and makes a GTF_DATA object that contains all the GTF file data.
 * This GTF_DATA is intended to be eventually indexed and used as input of one
 * of the library functions. This function is visible from outside the library.
 *
 * Parameters:
 * 		input:	a GTF file name (gzipped or not) or "-" for standard input
 *
 * Returns:		a GTF_DATA structure filled with GTF data
 */
__attribute__ ((visibility ("default")))
GTF_DATA *load_GTF(char *input) {
	/*
	 * a buffer to store rows of GTF file as they are read
	 */
	char *buffer = (char *)calloc(10000, sizeof(char));
	char *trim_buffer;

	/*
	 * string tables to split GTF rows and attributes into fields
	 */
	char **token, **attr;

	/*
	 * A loop integer
	 */
	int i;

	/*
	 * a counter for the rows
	 */
	int nb_row = 0;

	/*
	 * creates a GTF_READER to read from input
	 */
	TEXTFILE_READER *gr = get_gtf_reader(input);

	if (gr == NULL) return NULL;

	/*
	 * allocates memory for the GTF_DATA structure to return
	 */
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));

	/*
	 * allocates memory for one row in the table.
	 */
	ret->data = (GTF_ROW **)calloc(1, sizeof(GTF_ROW *));

	GTF_ROW *row, *previous_row;

	/*
	 * make the GTF column model
	 */
	make_columns();

	/*
	 * loop on the GTF file rows
	 */
	while (get_next_gtf_line(gr, buffer) != NULL) {
		if (*buffer != '#') {
			*(buffer + strlen(buffer) - 1) = 0;
			trim_buffer = trim_ip(buffer);

			/*
			 * allocates memory for a new row in the row list
			 */
			row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));
			if (nb_row == 0) ret->data[0] = row;

			/*
			 * split the row and check the number of fields (should be 9)
			 */
			i = split_ip(&token, trim_buffer, "\t");
			if (i != nb_column) {
				if (!strcmp(gr->filename, "-"))
					fprintf(stderr, "ERROR : standard input is not a valid GTF stream\n");
				else
					fprintf(stderr, "ERROR : GTF file %s is not valid (%d columns)\n", gr->filename, i);

				exit(0);
			}

			/*
			 * allocates memory for the 8 first fields
			 */
			row->field = (char **)calloc(8, sizeof(char *));

			/*
			 * store the 8 first fields
			 */
			for (i = 0; i < 8; i++) row->field[i] = strdup(token[i]);

			/*
			 * split the attributes
			 */
			row->attributes.nb = split_ip(&attr, token[8], ";");

			/*
			 * allocates memory for the attributes
			 */
			//row->attributes.attr = (ATTRIBUTE **)calloc(row->attributes.nb, sizeof(ATTRIBUTE *));
			row->attributes.attr = (ATTRIBUTE *)calloc(row->attributes.nb, sizeof(ATTRIBUTE));

			/*
			 * store the attributes
			 */
			for (i = 0; i < row->attributes.nb; i++) {
				//row->attributes.attr[i] = (ATTRIBUTE *)calloc(1, sizeof(ATTRIBUTE));
				//split_key_value(attr[i], &(row->attributes.attr[i]->key), &(row->attributes.attr[i]->value));
				split_key_value(attr[i], &((row->attributes.attr + i)->key), &((row->attributes.attr + i)->value));
			}

			/*
			 * keep the rank number of the row
			 */
			row->rank = nb_row;

			/*
			 * insert the new row in the list
			 */
			if (nb_row > 0) previous_row->next = row;
			previous_row = row;

			/*
			 * one more row has been read
			 */
			nb_row++;

			/*
			 * free memory used to split the row and the attributes
			 */
			free(token);
			free(attr);
		}
		else if (!strncmp(buffer, "##gff-version 3", 15)) {
			free(ret->data);
			free(ret);
			fprintf(stderr, "GFF3 format is not supported by libgtftk !\n");
			return NULL;
		}
	}

	/*
	 * store the final number of rows
	 */
	ret->size = nb_row;

	/*
	 * update the row table
	 */
	update_row_table(ret);

	/*
	 * free the buffer used to read GTF file
	 */
	free(buffer);
	if (gr != NULL) {
		free(gr->filename);
		free(gr);
	}
	return ret;
}

/*
 * This function returns the list of attributes in a GTF_DATA. It is used in
 * get_list function and in extract_data
 *
 * Parameters:
 * 		input:	the input GTF data
 *
 * Returns:		a STRING_LIST containing the list of attributes
 */
STRING_LIST *get_all_attributes(GTF_DATA *gtf_data) {
	/*
	 * Some loop variables
	 */
	int i, j;

	/*
	 * The list of attributes to return. Each time a new attribute is found,
	 * it will be added in this list.
	 */
	STRING_LIST *sl = (STRING_LIST *)calloc(1, sizeof(STRING_LIST));

	/*
	 * a variable to point on each attribute we want to look for in the
	 * result list
	 */
	char **pkey;

	/*
	 * A convenient variable to iterate on the rows of the input GTF data
	 */
	GTF_ROW *row;

	/*
	 * The loop on the rows
	 */
	for (j = 0; j < gtf_data->size; j++) {
		/*
		 * The current row
		 */
		row = gtf_data->data[j];

		/*
		 * The loop on the attributes of the current row
		 */
		for (i = 0; i < row->attributes.nb; i++) {

			/*
			 * the attribute to look for in the result list
			 */
			//pkey = &(row->attributes.attr[i]->key);
			pkey = &((row->attributes.attr + i)->key);

			/*
			 * If the attribute is not found in the result list, we add it with
			 * the lsearch function
			 */
			if (!lfind(pkey, sl->list, (size_t *)(&sl->nb), sizeof(char *), compatt)) {
				sl->list = (char **)realloc(sl->list, (sl->nb + 1) * sizeof(char *));
				lsearch(pkey, sl->list, (size_t *)(&sl->nb), sizeof(char *), compatt);
			}
		}
	}
	return sl;
}

/*
 * an empty function used to destroy the indexes when needed (tdestroy)
 */
void action_destroy(void *nodep) {

}

/*
 * this function is used to sort ROW_LIST elements of a new index because
 * the rows of the data files are randomized when loaded to avoid a linear
 * research tree
 */
static void action_sort(const void *nodep, const VISIT which, const int depth) {
	/*
	 * the row list of the current gene
	 */
	ROW_LIST *datap = *((ROW_LIST **)nodep);

	switch (which) {
			case preorder:
				break;

			/*
			 * The operations are made on internal nodes and leaves of the tree.
			 */
			case leaf :
			case postorder:
				qsort(datap->row, datap->nb_row, sizeof(int), comprow);
				break;

			case endorder:
				break;
		}
}

/*
 * This function is used to get an index on a column_name or an attribute.
 *
 * Parameters:
 * 		gtf_data:	the input GTF data
 * 		key:		a column name or an attribute name for wich the index is
 * 					looked for
 *
 * Returns:			an INDEX_ID structure containing the column number and the
 * 					index rank corresponding to to searched index.
 * 					If key is not a column name nor an attribute name, the
 * 					returned column number is set to -1.
 * 					If the index is not found, index rank is set to -1
 */
INDEX_ID *get_index(GTF_DATA *gtf_data, char *key) {
	int c, k;
	INDEX_ID *ret = calloc(1, sizeof(INDEX_ID));

	ret->column = ret->index_rank = -1;

	/*
	 * look for key in the 8 first column names and index the column if found
	 */
	for (c = 0; c < (nb_column - 1); c++)
		if (!strcmp(column[c]->name, key)) {
			/*
			 * key is the name of a column
			 * if an index on this column exists and correspond to the gtf_data
			 * parameter, we return it
			 */
			ret->column = c;
			for (k = 0; k < column[c]->nb_index; k++)
				if ((column[c]->index[k]->data != NULL) && (column[c]->index[k]->gtf_data == gtf_data)) {
					ret->index_rank = k;
					break;
				}
			break;
		}

	if (ret->column == -1) {
		/*
		 * key is not a column name so look for key in the attribute index
		 * table of the 8th column
		 */
		ret->column = 8;
		for (c = 0; c < column[8]->nb_index; c++)
			/*
			 * if an index on this parameter exists and correspond to the
			 * gtf_data parameter, we return it
			 */
			if (!strcmp(column[8]->index[c]->key, key)) {
				/*
				 * key is the name of a parameter
				 * if an index on this parameter exists and correspond to the
				 * gtf_data parameter, we return it
				 */
				if ((column[8]->index[c]->data != NULL) && (column[8]->index[c]->gtf_data == gtf_data)) {
					ret->index_rank = c;
					break;
				}
			}
	}
	return ret;
}

/*
 * This function check for an index on a GTF data with a column or a attribute.
 * If the index does not exist, it is created. Note that all the indexes are
 * stored in the column model (column external variable at the begining of this
 * file), declared in the libgtftk.c file.
 *
 * Parameters:
 * 		gtf_data:	the input GTF data
 * 		key:		the column or attribute name on which the index is to be made
 *
 * Returns:			the INDEX_ID of the requested index
 *
 */
INDEX_ID *index_gtf(GTF_DATA *gtf_data, char *key) {
	int j, k;

	srand(time(NULL));

	/*
	 * ask for the requested index
	 */
	INDEX_ID *index_id = get_index(gtf_data, key);

	/*
	 * The requested index does not exist
	 */
	if (index_id->index_rank == -1) {
		/*
		 * make a new index and add it in the corresponding column
		 */
		make_index(index_id, key);

		/*
		 * index is to be made on a column
		 */
		if (index_id->column != 8) {
			/*
			 * indexes all the rows of the GTF data with a column and put a
			 * reference on the GTF data in the INDEX
			 */
			for (k = 0; k < gtf_data->size; k++)
				index_row(k, gtf_data->data[k]->field[index_id->column], column[index_id->column]->index[index_id->index_rank]);
			column[index_id->column]->index[index_id->index_rank]->gtf_data = gtf_data;
		}
		/*
		 * index is to be made on an attribute
		 */
		else {
			/*
			 * indexes all the rows of the GTF data with an attribute and put a
			 * reference on the GTF data in the INDEX
			 */
			int *kk = (int *)calloc(gtf_data->size, sizeof(int)), kq, kr;
			for (k = 0; k < gtf_data->size; k++) kk[k] = k;
			for (k = 0; k < gtf_data->size; k++) {
				kr = k + rand() / (RAND_MAX / (gtf_data->size - k) + 1);
				kq = kk[k];
				kk[k] = kk[kr];
				kk[kr] = kq;
			}
			for (k = 0; k < gtf_data->size; k++) {
				for (j = 0; j < gtf_data->data[kk[k]]->attributes.nb; j++)
					/*if (!strcmp(key, gtf_data->data[kk[k]]->attributes.attr[j]->key)) {
						index_row(kk[k], gtf_data->data[kk[k]]->attributes.attr[j]->value, column[index_id->column]->index[index_id->index_rank]);
						break;
					}*/
					if (!strcmp(key, (gtf_data->data[kk[k]]->attributes.attr + j)->key)) {
						index_row(kk[k], (gtf_data->data[kk[k]]->attributes.attr + j)->value, column[index_id->column]->index[index_id->index_rank]);
						break;
					}
			}
			column[index_id->column]->index[index_id->index_rank]->gtf_data = gtf_data;
			twalk(column[index_id->column]->index[index_id->index_rank]->data, action_sort);
		}
	}
	return index_id;
}

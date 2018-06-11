/*
 * att_attributes.c
 *
 *  Created on: Jun 7, 2017
 *      Author: fafa
 */

#include "libgtftk.h"

/*
 * external functions declaration
 */
extern INDEX_ID *index_gtf(GTF_DATA *gtf_data, char *key);
extern int compare_row_list(const void *p1, const void *p2);
extern GTF_DATA *clone_gtf_data(GTF_DATA *gtf_data);
extern void add_attribute(GTF_ROW *row, char *key, char *value);

/*
 * global variables declaration
 */
extern COLUMN **column;

__attribute__ ((visibility ("default")))
GTF_DATA *add_attributes(GTF_DATA *gtf_data, char *features, char *key, char *new_key, char *inputfile_name) {
	int i;
	/*
	 * allocating memory for the GTF_DATA structure to return
	 */
	GTF_DATA *ret = clone_gtf_data(gtf_data);

	/*
	 * indexing the gtf with key
	 */
	INDEX_ID *ix = index_gtf(ret, key);

	GTF_ROW *row;

	FILE *input = fopen(inputfile_name, "ro");
	size_t buffersize = 1000;
	char *buffer = (char *)calloc(buffersize, sizeof(char));
	char *value, *new_value;
	ROW_LIST **find_row_list, *test_row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));

	while (getline(&buffer, &buffersize, input) > 0) {
		value = buffer;
		new_value = strchr(buffer, '\t') + 1;
		if (*(new_value + strlen(new_value) - 1) == '\n') *(new_value + strlen(new_value) - 1) = 0;
		*strchr(buffer, '\t') = 0;
		test_row_list->token = value;
		find_row_list = (ROW_LIST **)tfind(test_row_list, &(column[ix->column]->index[ix->index_rank]->data), compare_row_list);
		if (find_row_list != NULL)
			for (i = 0; i < (*find_row_list)->nb_row; i++) {
				row = ret->data[(*find_row_list)->row[i]];
				if (!strcmp(features, "*") || strstr(features, row->field[2]))
					add_attribute(row, new_key, new_value);
			}
	}
	if (test_row_list != NULL) {
		if (test_row_list->row != NULL) free(test_row_list->row);
		free(test_row_list);
	}
	free(buffer);
	fclose(input);
	return ret;
}

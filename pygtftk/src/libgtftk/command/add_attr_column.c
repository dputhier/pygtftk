/*
 * add_attr_column.c
 *
 *  Created on: Nov 14, 2017
 *      Author: dputhier (inspired from fafa code...)
 *      Purpose: Simply add a new column of attribute/value to each line.
 *      inputfile_name is a single column file with nrow(GTF_DATA) lines.
 *      File lines should be sorted according to the GTF_DATA object.
 *		Lines for which a new key should not be added should contain "?".
 */

#include "libgtftk.h"

/*
 * external functions declaration
 *
 */
extern GTF_DATA *clone_gtf_data(GTF_DATA *gtf_data);
extern void add_attribute(GTF_ROW *row, char *key, char *value);

/*
 * global variables declaration
 */

__attribute__ ((visibility ("default")))
GTF_DATA *add_attr_column(GTF_DATA *gtf_data, char *inputfile_name, char *new_key) {

	int i;
	GTF_DATA *ret = clone_gtf_data(gtf_data);
	GTF_ROW *row;
	FILE *input = fopen(inputfile_name, "ro");
	size_t buffersize = 1000;
	char *buffer = (char *)calloc(buffersize, sizeof(char));
	char *value;

	i = 0;
	while (getline(&buffer, &buffersize, input) > 0) {
		value = buffer;
		if (*(value + strlen(value) - 1) == '\n') *(value + strlen(value) - 1) = 0;
		row = ret->data[i];
		if (strcmp(value, "?")) add_attribute(row, new_key, value);
		i++;
	}
	free(buffer);
	return ret;
}

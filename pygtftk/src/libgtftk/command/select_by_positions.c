/*
 * select_by_pos.c
 *
 *  Created on: Mar 20, 2017
 *      Author: puthier (based on Fafa code...)
 *  Objective: select a set of line based on index/position
 */

#include "libgtftk.h"

extern int update_row_table(GTF_DATA *gtf_data);

__attribute__ ((visibility ("default")))
GTF_DATA *select_by_positions(GTF_DATA *gtf_data, int *pos, int size) {
	int i, j	;
	GTF_ROW *row, *previous_row = NULL;

	/*
	 * reserve memory for the GTF_DATA structure to return
	 */
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));

	// The size of the return gtf_data is the number of requested rows.
	ret->size = size;

	ret->data = (GTF_ROW **)calloc(1, sizeof(GTF_ROW *));

	for (i = 0; i < ret->size; i++) {
		row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));
		if (i == 0) ret->data[0] = row;
		row->rank = gtf_data->data[pos[i]]->rank;
		row->attributes.nb = gtf_data->data[pos[i]]->attributes.nb;
		row->field = (char **)calloc(8, sizeof(char*));
		for (j = 0; j < 8; j++) row->field[j] = strdup(gtf_data->data[pos[i]]->field[j]);
		row->attributes.nb = gtf_data->data[pos[i]]->attributes.nb;
		/*row->attributes.attr = (ATTRIBUTE **)calloc(row->attributes.nb, sizeof(ATTRIBUTE *));
		for (j = 0; j < row->attributes.nb; j++) {
			row->attributes.attr[j] = (ATTRIBUTE *)calloc(1, sizeof(ATTRIBUTE));
			row->attributes.attr[j]->key = strdup(gtf_data->data[pos[i]]->attributes.attr[j]->key);
			row->attributes.attr[j]->value = strdup(gtf_data->data[pos[i]]->attributes.attr[j]->value);
			if (j > 0) row->attributes.attr[j - 1]->next = row->attributes.attr[j];
		}*/
		row->attributes.attr = (ATTRIBUTE *)calloc(row->attributes.nb, sizeof(ATTRIBUTE));
		for (j = 0; j < row->attributes.nb; j++) {
			(row->attributes.attr + j)->key = strdup((gtf_data->data[pos[i]]->attributes.attr + j)->key);
			(row->attributes.attr + j)->value = strdup((gtf_data->data[pos[i]]->attributes.attr + j)->value);
		}
		if (i > 0) previous_row->next = row;
		previous_row = row;
	}

	update_row_table(ret);

	return ret;
}

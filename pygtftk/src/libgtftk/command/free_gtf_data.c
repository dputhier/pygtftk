/*
 * free_gtf_data.c
 *
 *  Created on: Mar 20, 2017
 *      Author: fafa
 */

#include "libgtftk.h"

/*
 * global variables in libgtftk.c
 */
extern COLUMN **column;
extern int nb_column;

int update_index_table(COLUMN *col) {
	INDEX *idx;
	int i;

	if (col->index != NULL) {
		idx = col->index[0];
		col->index = (INDEX **)realloc(col->index, col->nb_index * sizeof(INDEX *));
		for (i = 0; i < col->nb_index; i++) {
			col->index[i] = idx;
			idx = idx->next;
		}
	}
	return 0;
}

/*
 * This function is called by twalk on an index of ROW_LIST to destroy all the
 * elements of the tree.
 */
static void destroy_row_list_tree(const void *nodep, const VISIT which, const int depth) {
	ROW_LIST *rl = (ROW_LIST *)nodep;

	switch (which) {
		case preorder:
			break;
		case postorder:
			break;
		case endorder:
		case leaf:
			if (rl->token != NULL) free(rl->token);
			if (rl->row != NULL) free(rl->row);
			break;
	}
}

__attribute__ ((visibility ("default")))
int free_gtf_data(GTF_DATA *gtf_data) {
	int i, j, c;
	GTF_ROW *row;
	INDEX *pindex, *pindex0;
	ROW_LIST *row_list;

	row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));
	row_list->token = strdup("*");

	if (gtf_data != NULL) {
		for (i = 0; i < gtf_data->size; i++) {
			row = gtf_data->data[i];
			for (j = 0; j < 8; j++) if (row->field[j] != NULL) free(row->field[j]);
			free(row->field);

			for (j = 0; j < row->attributes.nb; j++) {
				//if (row->attributes.attr[j]->key != NULL) free(row->attributes.attr[j]->key);
				//if (row->attributes.attr[j]->value != NULL) free(row->attributes.attr[j]->value);
				//free(row->attributes.attr[j]);
				if ((row->attributes.attr + j)->key != NULL) free((row->attributes.attr + j)->key);
				if ((row->attributes.attr + j)->value != NULL) free((row->attributes.attr + j)->value);
			}
			free(row->attributes.attr);
			free(row);
		}
		free(gtf_data->data);
		gtf_data->data = NULL;
		for (c = 0; c < nb_column; c++) {
			if (column[c]->nb_index != 0)
				pindex = column[c]->index[0];
			else
				pindex = NULL;
			pindex0 = NULL;
			while (pindex != NULL) {
				if (pindex->gtf_data == gtf_data) {
					twalk(pindex->data, destroy_row_list_tree);
					free(pindex->key);
					column[c]->nb_index--;
					if (pindex0 == NULL) {
						pindex0 = pindex->next;
						free(pindex);
						if (pindex == column[c]->index[0]) column[c]->index[0] = pindex0;
						pindex = pindex0;
						pindex0 = NULL;
					}
					else {
						pindex0->next = pindex->next;
						free(pindex);
						if (pindex == column[c]->index[0]) column[c]->index[0] = pindex0->next;
						pindex = pindex0->next;
					}
					update_index_table(column[c]);
				}
				else {
					pindex0 = pindex;
					pindex = pindex->next;
				}
			}
			update_index_table(column[c]);
		}
		free(gtf_data);
		gtf_data = NULL;
	}
	free(row_list->token);
	free(row_list);
	return 0;
}

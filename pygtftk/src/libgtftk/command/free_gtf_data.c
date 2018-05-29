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

extern void *bookmem(int nb, int size, char *file, const char *func, int line);
extern void *rebookmem(void *ptr, int size, char *file, const char *func, int line);
extern void freemem(void *ptr, char *file, const char *func, int line);
extern char *dupstring(const char *s, char *file, const char *func, int line);

int update_index_table(COLUMN *col) {
	INDEX *idx;
	int i;

	if (col->index != NULL) {
		idx = col->index[0];
		//col->index = (INDEX **)realloc(col->index, col->nb_index * sizeof(INDEX *));
		col->index = (INDEX **)rebookmem(col->index, col->nb_index * sizeof(INDEX *), __FILE__, __func__, __LINE__);
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
			if (rl->token != NULL) freemem(rl->token, __FILE__, __func__, __LINE__);
			if (rl->row != NULL) freemem(rl->row, __FILE__, __func__, __LINE__);
			break;
	}
}

__attribute__ ((visibility ("default")))
int free_gtf_data(GTF_DATA *gtf_data) {
	int i, j, c;
	GTF_ROW *row;
	INDEX *pindex, *pindex0;
	ROW_LIST *row_list;

	//row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));
	row_list = (ROW_LIST *)bookmem(1, sizeof(ROW_LIST), __FILE__, __func__, __LINE__);
	//row_list->token = strdup("*");
	row_list->token = dupstring("*", __FILE__, __func__, __LINE__);

	if (gtf_data != NULL) {
		for (i = 0; i < gtf_data->size; i++) {
			row = gtf_data->data[i];
			for (j = 0; j < 8; j++) if (row->field[j] != NULL) freemem(row->field[j], __FILE__, __func__, __LINE__);
			freemem(row->field, __FILE__, __func__, __LINE__);

			for (j = 0; j < row->attributes.nb; j++) {
				if (row->attributes.attr[j]->key != NULL) freemem(row->attributes.attr[j]->key, __FILE__, __func__, __LINE__);
				if (row->attributes.attr[j]->value != NULL) freemem(row->attributes.attr[j]->value, __FILE__, __func__, __LINE__);
				freemem(row->attributes.attr[j], __FILE__, __func__, __LINE__);
			}
			freemem(row->attributes.attr, __FILE__, __func__, __LINE__);
			freemem(row, __FILE__, __func__, __LINE__);
		}
		freemem(gtf_data->data, __FILE__, __func__, __LINE__);
		gtf_data->data = NULL;
		//fprintf(stderr, "freed %d\n", gtf_data->size);
		for (c = 0; c < nb_column; c++) {
			//fprintf(stderr, "  col = %s %d\n", column[c]->name, column[c]->nb_index);
			if (column[c]->nb_index != 0)
				pindex = column[c]->index[0];
			else
				pindex = NULL;
			pindex0 = NULL;
			while (pindex != NULL) {
				if (pindex->gtf_data == gtf_data) {
					//fprintf(stderr, "    freeing index %s ...", pindex->key);
					twalk(pindex->data, destroy_row_list_tree);
					//fprintf(stderr, " OK\n");
					freemem(pindex->key, __FILE__, __func__, __LINE__);
					column[c]->nb_index--;
					if (pindex0 == NULL) {
						pindex0 = pindex->next;
						freemem(pindex, __FILE__, __func__, __LINE__);
						if (pindex == column[c]->index[0]) column[c]->index[0] = pindex0;
						pindex = pindex0;
						pindex0 = NULL;
					}
					else {
						pindex0->next = pindex->next;
						freemem(pindex, __FILE__, __func__, __LINE__);
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
			//fprintf(stderr, "Updating index table ...");
			update_index_table(column[c]);
			//fprintf(stderr, " OK\n");
		}
		freemem(gtf_data, __FILE__, __func__, __LINE__);
		gtf_data = NULL;
	}
	freemem(row_list->token, __FILE__, __func__, __LINE__);
	freemem(row_list, __FILE__, __func__, __LINE__);
	return 0;
}

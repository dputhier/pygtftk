/*
 * get_list.c
 *
 *  Created on: Apr 14, 2017
 *      Author: fafa
 */

#include "libgtftk.h"

extern INDEX_ID *index_gtf(GTF_DATA *gtf_data, char *key);
extern COLUMN **column;
extern STRING_LIST *get_all_attributes(GTF_DATA *gtf_data);
extern TTEXT *vret;

static void action_list(const void *nodep, const VISIT which, const int depth) {
	ROW_LIST *datap;
	char tmp[100];

	switch (which) {
		case preorder:
			break;
		case postorder:
		case leaf:
			datap = *((ROW_LIST **)nodep);
			vret->data = (char ***)realloc(vret->data, (vret->size + 1) * sizeof(char **));
			vret->data[vret->size] = (char **)calloc(2, sizeof(char *));
			sprintf(tmp, "%d", datap->nb_row);
			vret->data[vret->size][0] = strdup(tmp);
			vret->data[vret->size][1] = strdup(datap->token);
			vret->size++;
			break;
		case endorder:
			break;
	}
}

__attribute__ ((visibility ("default")))
TTEXT *get_feature_list(GTF_DATA *gtf_data) {
	/*
	 * reserve memory for the TTEXT structure to return
	 */
	vret = (TTEXT *)calloc(1, sizeof(TTEXT));

	/*
	 * indexing the GTF_DATA with the feature column
	 */
	INDEX_ID *index_id = index_gtf(gtf_data, "feature");

	/*
	 * tree browsing of the feature index
	 */
	twalk(column[index_id->column]->index[index_id->index_rank]->data, action_list);

	return vret;
}

__attribute__ ((visibility ("default")))
TTEXT *get_seqid_list(GTF_DATA *gtf_data) {
	/*
	 * reserve memory for the TTEXT structure to return
	 */
	vret = (TTEXT *)calloc(1, sizeof(TTEXT));

	/*
	 * indexing the GTF_DATA with the seqid column
	 */
	INDEX_ID *index_id = index_gtf(gtf_data, "seqid");

	/*
	 * tree browsing of the feature index
	 */
	twalk(column[index_id->column]->index[index_id->index_rank]->data, action_list);

	return vret;
}

__attribute__ ((visibility ("default")))
TTEXT *get_attribute_list(GTF_DATA *gtf_data) {
	int i;
	STRING_LIST *sl = get_all_attributes(gtf_data);

	/*
	 * reserve memory for the TTEXT structure to return
	 */
	vret = (TTEXT *)calloc(1, sizeof(TTEXT));

	/*
	 * reserve memory for the data table (a list in this case)
	 */
	vret->data = (char ***)calloc(sl->nb, sizeof(char **));

	/*
	 * fill the list-like table
	 */
	for (i = 0; i < sl->nb; i++) {
		vret->data[i] = (char **)calloc(1, sizeof(char *));
		vret->data[i][0] = strdup(sl->list[i]);
	}
	/*
	 * setup the number of strings
	 */
	vret->size = sl->nb;

	return vret;
}

__attribute__ ((visibility ("default")))
TTEXT *get_attribute_values_list(GTF_DATA *gtf_data, char *attribute) {
	/*
	 * reserve memory for the TTEXT structure to return
	 */
	vret = (TTEXT *)calloc(1, sizeof(TTEXT));

	/*
	 * indexing the GTF_DATA with the provided attribute
	 */
	INDEX_ID *index_id = index_gtf(gtf_data, attribute);

	/*
	 * tree browsing of the feature index
	 */
	twalk(column[index_id->column]->index[index_id->index_rank]->data, action_list);

	return vret;
}

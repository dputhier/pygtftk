/*
 * select_by_transcript_size.c
 *
 *  Created on: Jan 12, 2017
 *      Author: fafa
 *
 * Implementation of select_by_transcript_size function.
 * This function select rows in GTF_DATA that correspond to transcripts with a
 * size between the given min and max values. The size of a transcript is
 * the sum of the sizes of his exons.
 */

#include "libgtftk.h"

/*
 * external functions declaration
 */
extern INDEX_ID *index_gtf(GTF_DATA *gtf_data, char *key);
extern int comprow(const void *m1, const void *m2);
extern int add_row_list(ROW_LIST *src, ROW_LIST *dst);
extern int update_row_table(GTF_DATA *gtf_data);
extern void add_attribute(GTF_ROW *row, char *key, char *value);

/*
 * global variables declaration
 */
extern COLUMN **column;
extern ROW_LIST *row_list;
extern GTF_DATA *gtf_d;
extern int min_ts, max_ts;

/*
 * The function used by twalk on each node of the index tree.
 * This function is used to browse an index on "transcript_id" attribute
 * containing ROW_LIST elements. For each transcript, we make the sum of his
 * exon lengthes and if this sum is between min and max values, the associated
 * ROW_LIST element is merged with the local row_list.
 * For information about the parameters, see man pages of twalk.
 */
static void action_sbts(const void *nodep, const VISIT which, const int depth) {
	ROW_LIST *datap;
	GTF_ROW *row;
	int i, trsize;

	switch (which) {
		case preorder:
			break;

		/*
		 * The operations are made on internal nodes and leaves of the tree.
		 */
		case postorder:
		case leaf:
			datap = *((ROW_LIST **)nodep);

			// computing the size of the transcript (the sum of exon sizes)
			trsize = 0;
			for (i = 0; i < datap->nb_row; i++) {
				row = gtf_d->data[datap->row[i]];
				if (!strcmp(row->field[2], "exon"))
					trsize += (atoi(row->field[4]) - atoi(row->field[3]) + 1);
			}

			/* if this size is between min and max, we add the transcript rows
			 * in row_list
			 */
			if ((trsize >= min_ts) && (trsize <= max_ts)) add_row_list(datap, row_list);
			break;

		case endorder:
			break;
	}
}

/*
 * select_by_transcript_size function selects rows in GTF_DATA that correspond
 * to transcripts with a size between min and max values.
 *
 * Parameters:
 * 		gtf_data:	a GTF_DATA structure
 * 		min:		the minimum size of transcripts
 * 		max:		the maximum size of transcripts
 *
 * Returns:			a GTF_DATA structure that contains the result of the query
 */
__attribute__ ((visibility ("default")))
GTF_DATA *select_by_transcript_size(GTF_DATA *gtf_data, int min, int max) {
	int i, k;
	INDEX_ID *trid_index_id;

	/*
	 * reserve memory for the GTF_DATA structure to return
	 */
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));

	/*
	 * reserve memory for the local ROW_LIST
	 */
	row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));

	/*
	 * setup local variables to allow action_sbts function to access these
	 * values
	 */
	gtf_d = gtf_data;
	min_ts = min;
	max_ts = max;

	/*
	 * indexes the GTF_DATA with transcript_id attribute to create an index
	 * containing, for each transcript, the list of his related rows. The real
	 * rank of the index (- 8) in the attributes column index table is stored
	 * in i.
	 */
	trid_index_id = index_gtf(gtf_data, "transcript_id");

	// tree browsing of the transcript_id index
	twalk(column[trid_index_id->column]->index[trid_index_id->index_rank]->data, action_sbts);

	/*
	 * we sort the resulting row list to be sure to respect the original order
	 * of the transcripts
	 */
	qsort(row_list->row, row_list->nb_row, sizeof(int), comprow);

	/*
	 * now we fill the resulting GTF_DATA with the found rows and return it
	 */
	ret->data = (GTF_ROW **)calloc(row_list->nb_row, sizeof(GTF_ROW *));
	GTF_ROW *row, *previous_row = NULL;
	for (i = 0; i < row_list->nb_row; i++) {
		row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));
		row->field = (char **)calloc(8, sizeof(char *));
		if (i == 0) ret->data[0] = row;
		for (k = 0; k < gtf_data->data[row_list->row[i]]->attributes.nb; k++)
			//add_attribute(row, gtf_data->data[row_list->row[i]]->attributes.attr[k]->key,
			//				gtf_data->data[row_list->row[i]]->attributes.attr[k]->value);
			add_attribute(row, (gtf_data->data[row_list->row[i]]->attributes.attr + k)->key,
							(gtf_data->data[row_list->row[i]]->attributes.attr + k)->value);
		for (k = 0; k < 8; k++)
			row->field[k] = strdup(gtf_data->data[row_list->row[i]]->field[k]);
		row->rank = gtf_data->data[row_list->row[i]]->rank;
		if (i > 0) previous_row->next = row;
		previous_row = row;
	}
	ret->size = row_list->nb_row;
	update_row_table(ret);
	return ret;
}

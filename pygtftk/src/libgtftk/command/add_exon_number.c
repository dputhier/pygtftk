/*
 * add_exon_number.c
 *
 *  Created on: Jul 31, 2017
 *      Author: fafa
 */

#include "libgtftk.h"

/*
 * external functions declaration
 */
extern void add_attribute(GTF_ROW *row, char *key, char *value);
extern INDEX_ID *index_gtf(GTF_DATA *gtf_data, char *key);
extern GTF_DATA *clone_gtf_data(GTF_DATA *gtf_data);

/*
 * global variables declaration
 */
extern COLUMN **column;
extern GTF_DATA *gtf_d;
extern SORT_ROW *sort_row;
extern int nb_sort_row;
extern char *enf;

/*
 * the comparison function for SORT_ROW structures. Used to sort the exon rows
 * of each transcript to get their rank, based on the genomic locations
 */
int compare_sort_row(const void *r1, const void *r2) {
	SORT_ROW *sr1 = (SORT_ROW *)r1;
	SORT_ROW *sr2 = (SORT_ROW *)r2;
	return sr1->value - sr2->value;
}

/*
 * The tree parsing function used by twalk.
 * This function is used to browse an index on "transcript_id" attribute
 * containing ROW_LIST elements. For each exon line, we add an attribute called
 * exon_number to indicate the exon's rank in the transcript.
 */
static void action_aen(const void *nodep, const VISIT which, const int depth) {
	ROW_LIST *datap;
	GTF_ROW *row;
	char tmp[10];
	int i, start, end;

	switch (which) {
		case preorder:
			break;

		case leaf:
		case postorder:
			datap = *((ROW_LIST **)nodep);
			nb_sort_row = 0;
			for (i = 0; i < datap->nb_row; i++) {
				row = gtf_d->data[datap->row[i]];
				if (!strcmp(row->field[2], "exon")) {
					nb_sort_row++;
					sort_row = (SORT_ROW *)realloc(sort_row, nb_sort_row * sizeof(SORT_ROW));
					start = atoi(row->field[3]);
					end = atoi(row->field[4]);
					sort_row[nb_sort_row - 1].row = i;
					if (*(row->field[6]) == '+')
						sort_row[nb_sort_row - 1].value = start;
					else
						sort_row[nb_sort_row - 1].value = -end;
				}
			}
			qsort(sort_row, nb_sort_row, sizeof(SORT_ROW), compare_sort_row);
			for (i = 0; i < nb_sort_row; i++) {
				row = gtf_d->data[datap->row[sort_row[i].row]];
				sprintf(tmp, "%d", i + 1);
				add_attribute(row, enf, tmp);
			}
			break;

		case endorder:
			break;
	}
}

/*
 *
 */
__attribute__ ((visibility ("default")))
GTF_DATA *add_exon_number(GTF_DATA *gtf_data, char *exon_number_field) {
	/*
	 * copy the initial GTF data
	 */
	GTF_DATA *ret = clone_gtf_data(gtf_data);

	/*
	 * indexing the gtf with transcript_id
	 */
	INDEX_ID *ix = index_gtf(ret, "transcript_id");

	/*
	 * setup local variables to allow action_sbnoe function to access these
	 * values
	 */
	gtf_d = ret;
	nb_sort_row = 0;
	sort_row = NULL;
	if (exon_number_field != NULL)
		enf = strdup(exon_number_field);
	else
		enf = strdup("exon_number");

	/*
	 * tree browsing of the transcript_id index
	 */
	twalk(column[ix->column]->index[ix->index_rank]->data, action_aen);

	free(enf);

	return ret;
}

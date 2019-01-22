#include "libgtftk.h"

extern INDEX_ID *index_gtf(GTF_DATA *gtf_data, char *key);
extern int compare_row_list(const void *p1, const void *p2);
extern int comprow(const void *m1, const void *m2);
extern int add_row_list(ROW_LIST *src, ROW_LIST *dst);
extern int add_row(int row, ROW_LIST *dst);
extern int update_row_table(GTF_DATA *gtf_data);
extern void add_attribute(GTF_ROW *row, char *key, char *value);

/*
 * global variables declaration
 */
extern COLUMN **column;
extern INDEX_ID *tid_index;
extern ROW_LIST *row_list, *test_row_list, **find_row_list;
extern GTF_DATA *gtf_d;
extern int tr_type;

int string_cmp (const void *v1, const void *v2) {
	char *s1 = *(char **)v1;
	char *s2 = *(char **)v2;
	return strcmp(s1, s2);
}
int get_trid_list(ROW_LIST *rl, char ***tr_list) {
	int i, j;
	size_t n = 0;
	GTF_ROW *row;

	for (i = 0; i < rl->nb_row; i++) {
		row = gtf_d->data[rl->row[i]];
		for (j = 0; j < row->attributes.nb; j++)
			//if (!strcmp(row->attributes.attr[j]->key, "transcript_id")) {
			if (!strcmp((row->attributes.attr + j)->key, "transcript_id")) {
				*tr_list = (char **)realloc(*tr_list, (n + 1) * sizeof(char *));
				//lsearch(&(row->attributes.attr[j]->value), *tr_list, &n, sizeof(char *), string_cmp);
				lsearch(&((row->attributes.attr + j)->value), *tr_list, &n, sizeof(char *), string_cmp);
				break;
			}
	}

	return n;
}

static void action_st(const void *nodep, const VISIT which, const int depth) {
	// loop variables
	int i, k;

	// the row list of the current gene
	ROW_LIST *datap = *((ROW_LIST **)nodep);

	// the current row
	GTF_ROW *row;

	char **tr_list = NULL;
	int nb_tr;

	/*
	 * the feature value of a row
	 */
	char *feature;

	// the "gene" row rank
	int rg = -1;

	// the size of the current transcript
	int trsize;

	// the coordinates of the current exon or transcript
	int start, end;

	int min_i = 0, min_trsize;
	int max_i = 0, max_trsize;
	int most_5p, most_5p_i = 0;

	switch (which) {
		case preorder:
			break;
		case leaf :
		case postorder:
			max_trsize = 0;
			min_trsize = 10000000;

			/*
			 * loop on the rows of a gene
			 */
			for (i = 0; i < datap->nb_row; i++) {
				// the current row
				row = gtf_d->data[datap->row[i]];
				feature = row->field[2];
				if (!strcmp(feature, "gene"))
					rg = datap->row[i];
			}
			nb_tr = get_trid_list(datap, &tr_list);

			/*
			 * we set most_5p to a great number if the gene is on the positive strand
			 */
			most_5p = 0;
			if (datap->nb_row > 0) {
				row = gtf_d->data[rg];
				if (*(row->field[6]) == '+') most_5p = 300000000;
			}

			//fprintf(stderr, "nb_tr = %d\n", nb_tr);
			for (i = 0; i < nb_tr; i++) {
				row_list->token = tr_list[i];
				//fprintf(stderr, " %s\n", row_list->token);
				find_row_list = (ROW_LIST **)tfind(row_list, &(column[8]->index[tid_index->index_rank]->data), compare_row_list);
				if (find_row_list != NULL) {
					trsize = 0;
					/*
					 * we sort the resulting row list to respect the original order of the
					 * transcripts
					 */
					qsort((*find_row_list)->row, (*find_row_list)->nb_row, sizeof(int), comprow);

					/*
					 * we loop on all exon rows of the transcript
					 * to compute his length
					 */
					for (k = 0; k < (*find_row_list)->nb_row; k++) {
						row = gtf_d->data[(*find_row_list)->row[k]];
						if (!strcmp(row->field[2], "exon")) {
							start = atoi(row->field[3]);
							end = atoi(row->field[4]);
							trsize += (end - start + 1);
							if (*(row->field[6]) == '+') {
								if (most_5p > start) {
									most_5p = start;
									most_5p_i = i;
								}
							}
							else {
								if (most_5p < end) {
									most_5p = end;
									most_5p_i = i;
								}
							}
						}
					}

					/*
					 * now we check if it is the current shortest
					 * transcript and, if it is the case, we storelsearch(&(row->attributes.attr[j]->value), *tr_list, &n, sizeof(char *), string_cmp);:
					 * 	the length
					 * 	the rank of the first row
					 * 	the rank of the transcript_id attribute
					 */
					if (trsize < min_trsize) {
						min_trsize = trsize;
						min_i = i;
					}
					/*
					 * the same for the longest transcript
					 */
					if (trsize > max_trsize) {
						max_trsize = trsize;
						max_i = i;
					}
				}
			}

			/*
			 * after parsing the gene, min_i and max_i are the ranks of the
			 * first rows ("transcript" rows) of the shortest and longest
			 * transcripts of the gene and min_j and max_j are the ranks of the
			 * "transcript_id" attribute of min_i and max_i rows
			 */

			if (tr_type == SHORTEST_TRANSCRIPT) {
				/*
				 * we look for the shortest transcript
				 */
				test_row_list->token = tr_list[min_i];
			}
			else if (tr_type == LONGEST_TRANSCRIPT) {
				/*
				 * we look for the longest transcript
				 */
				test_row_list->token = tr_list[max_i];
			}
			else if (tr_type == MOST5P_TRANSCRIPT) {
				/*
				 * we look for the longest transcript
				 */
				test_row_list->token = tr_list[most_5p_i];
			}
			/*
			 * get transcript rows
			 */
			find_row_list = tfind(test_row_list, &(column[8]->index[tid_index->index_rank]->data), compare_row_list);

			/*
			 * we add the gene and transcript rows in the row_list structure
			 */
			add_row_list(*find_row_list, row_list);
			if (rg != -1) add_row(rg, row_list);
			if (nb_tr > 0) free(tr_list);
			break;

		case endorder:
			break;
	}
}

__attribute__ ((visibility ("default")))
GTF_DATA *select_transcript(GTF_DATA *gtf_data, int type) {
	INDEX_ID *gid_index;
	int i, k;
	tr_type = type;

	/*
	 * reserve memory for the GTF_DATA structure to return
	 */
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));

	/*
	 * indexes the GTF_DATA with gene_id and transcript_id attributes. The
	 * first index contains, for each gene, the list of his related rows.
	 * The rank of this index is 0. The second index contains, for each
	 * transcript, the list of his related rows. The rank of this index is 1.
	 */
	gid_index = index_gtf(gtf_data, "gene_id");
	tid_index = index_gtf(gtf_data, "transcript_id");

	/*
	 * setup local variables to allow action_sst function to access these
	 * values
	 */
	gtf_d = gtf_data;

	/*
	 * reserve memory for the local ROW_LIST
	 */
	row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));
	test_row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));

	// tree browsing of the gene_id index (rank 0)
	twalk(column[8]->index[gid_index->index_rank]->data, action_st);

	/*
	 * we sort the resulting row list to respect the original order of the
	 * transcripts
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
			//	gtf_data->data[row_list->row[i]]->attributes.attr[k]->value);
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

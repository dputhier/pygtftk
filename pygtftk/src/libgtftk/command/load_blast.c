/*
 * load_blast.c
 *
 *  Created on: Jul 17, 2018
 *      Author: fafa
 */

#include "libgtftk.h"

/*
 * external functions in blast_reader.c
 */
extern TEXTFILE_READER *get_blast_reader(char *query);
extern char *get_blast_header(TEXTFILE_READER *, BLAST_HEADER *);
extern char *get_next_blast_hsp(TEXTFILE_READER *, BLAST_HSP *, char *);

/*
 * external functions in column.c
 */
extern void make_columns(void);

/*
 * external functions in load_gtf.c
 */
extern int update_row_table(GTF_DATA *gtf_data);

/*
 * Create an attribute and add it in the GTF_ROW given in parameters
 */
void new_attribute(char *key, char *value, GTF_ROW *row) {
	row->attributes.nb++;
	//row->attributes.attr = (ATTRIBUTE **)realloc(row->attributes.attr, row->attributes.nb * sizeof(ATTRIBUTE *));
	//row->attributes.attr[row->attributes.nb - 1] = (ATTRIBUTE *)calloc(1, sizeof(ATTRIBUTE));
	//row->attributes.attr[row->attributes.nb - 1]->key = strdup(key);
	//row->attributes.attr[row->attributes.nb - 1]->value = strdup(value);
	row->attributes.attr = (ATTRIBUTE *)realloc(row->attributes.attr, row->attributes.nb * sizeof(ATTRIBUTE));
	(row->attributes.attr + row->attributes.nb - 1)->key = strdup(key);
	(row->attributes.attr + row->attributes.nb - 1)->value = strdup(value);
}

/*
 * Make a GTF_ROW from a BLASTN HSP
 */
GTF_ROW *make_row(BLAST_HSP *hsp, GTF_DATA *gtf_data, int rank) {
	char *buffer = (char *)calloc(10000, sizeof(char));

	GTF_ROW *row = (GTF_ROW *)calloc(1, sizeof(GTF_ROW));
	if (rank == 0) gtf_data->data[0] = row;
	row->field = (char **)calloc(8, sizeof(char *));
	row->field[0] = strdup(hsp->bs.subject_name);
	row->field[1] = strdup(hsp->bh.program_name);
	row->field[2] = strdup("HSP");
	sprintf(buffer, "%d", hsp->subject_start);
	row->field[3] = strdup(buffer);
	sprintf(buffer, "%d", hsp->subject_end);
	row->field[4] = strdup(buffer);
	sprintf(buffer, "%f", hsp->score);
	row->field[5] = strdup(buffer);
	row->field[6] = (char *)calloc(2, sizeof(char));
	*(row->field[6]) = hsp->strand_subject;
	row->field[7] = (char *)calloc(2, sizeof(char));
	*(row->field[7]) = '.';
	row->attributes.nb = 0;
	row->attributes.attr = NULL;
	new_attribute("database_name", hsp->bh.database_name, row);
	sprintf(buffer, "%u", hsp->bh.database_length);
	new_attribute("database_length", buffer, row);
	sprintf(buffer, "%d", hsp->bh.database_nb_sequences);
	new_attribute("database_nb_sequences", buffer, row);
	new_attribute("query_name", hsp->bq.query_name, row);
	sprintf(buffer, "%d", hsp->bq.query_length);
	new_attribute("query_length", buffer, row);
	sprintf(buffer, "%d", hsp->bs.subject_length);
	new_attribute("subject_length", buffer, row);
	sprintf(buffer, "%g", hsp->expect);
	new_attribute("expect", buffer, row);
	new_attribute("identities", hsp->identities, row);
	sprintf(buffer, "%d", hsp->identities_percent);
	new_attribute("identities_percent", buffer, row);
	if (hsp->gaps != NULL) {
		new_attribute("gaps", hsp->gaps, row);
		sprintf(buffer, "%d", hsp->gap_percent);
		new_attribute("gaps_percent", buffer, row);
	}
	row->rank = rank;

	free(buffer);
	return row;
}

void free_hsp(BLAST_HSP *hsp) {
	free(hsp->gaps);
	free(hsp->identities);
	free(hsp->bq.query_name);
	free(hsp->bs.subject_name);
	free(hsp->bh.database_name);
	free(hsp->bh.program_name);
}

/*
 * The load_blast function makes a GTF_DATA structure from a standard pairwise
 * BLASTN output file. The feature type is "HSP". Note that there is no gene_id
 * or transcript_id attributes in the results, so the GTF_DATA is not usable
 * with functions that need these attributes.
 */
__attribute__ ((visibility ("default")))
GTF_DATA *load_blast(char *input) {
	/*
	 * A reader for the BLAST file.
	 */
	TEXTFILE_READER *gr = get_blast_reader(input);

	/*
	 * If input is not a valid file or standard input, return NULL pointer
	 */
	if (gr == NULL) return NULL;

	/*
	 * Allocates memory for the returned GTF_DATA
	 */
	GTF_DATA *ret = (GTF_DATA *)calloc(1, sizeof(GTF_DATA));
	ret->data = (GTF_ROW **)calloc(1, sizeof(GTF_ROW *));

	/*
	 * Make the column model if not already done
	 */
	make_columns();

	/*
	 * A structure to read the HSPs is the BLAST file, one by one.
	 */
	BLAST_HSP *hsp = (BLAST_HSP *)calloc(1, sizeof(BLAST_HSP));

	/*
	 * Read the BLAST header of the file, i.e. the rows before the first HSP.
	 * The first row of the first HSP is returned to avoid going back in the
	 * BLAST file.
	 */
	char *r = get_blast_header(gr, &(hsp->bh));

	/*
	 * Some convenient local variables used to build the GTF rows.
	 */
	GTF_ROW *row, *previous_row;

	/*
	 * The row number in the future GTF_DATA
	 */
	int n = 0;

	/*
	 * The main loop on HSPs. The function get_next_blast_hsp returns the
	 * string just after the HSP that has been read. Then, we don't need to go
	 * back in the BLAST file to read the next HSP.
	 */
	while ((r = get_next_blast_hsp(gr, hsp, r)) != NULL) {

		/*
		 * Make a GTF_ROW from the HSP
		 */
		row = make_row(hsp, ret, n);

		/*
		 * Link the new row to the last one if needed
		 */
		if (n > 0) previous_row->next = row;
		previous_row = row;

		/*
		 * increment the row number
		 */
		n++;
	}

	/*
	 * Do the job for the last HSP
	 */
	row = make_row(hsp, ret, n);
	if (n > 0) previous_row->next = row;
	previous_row = row;
	n++;

	/*
	 * free memory in hsp
	 */
	free_hsp(hsp);

	/*
	 * Free some local memory allocations
	 */
	free(r);
	free(gr->filename);
	if (gr->plainfile != NULL) free(gr->plainfile);
	if (gr->gzfile != NULL) gzclose(gr->gzfile);
	free(gr);
	free(hsp);

	/*
	 * Update the size of the resuls
	 */
	ret->size = n;

	/*
	 * Make the table of HSPs from the linked list built above
	 */
	update_row_table(ret);

	return ret;
}

/*
 * blast_reader.c
 *
 *  Created on: Jul 17, 2018
 *      Author: fafa
 */

#include "libgtftk.h"

extern int split_ip(char ***tab, char *s, char *delim);
extern char *trim_ip(char *);

/*
 * Build and return a reader on a BLASTN output file.
 */
TEXTFILE_READER *get_blast_reader(char *query) {
	TEXTFILE_READER *gr = (TEXTFILE_READER *)calloc(1, sizeof(TEXTFILE_READER));
	char *tmp = (char *)calloc(1000, sizeof(char)), *query_filename = NULL;

	if (!access(query, F_OK) || !strcmp(query, "-"))
		query_filename = strdup(query);

	if (query_filename != NULL) {
		/*
		 * here we got a valid file or standard input in query_filename
		 */
		if (strstr(query_filename, ".gz")) {
			gr->filename = query_filename;
			gr->gz = 1;
			gr->gzfile = gzopen(gr->filename, "r");
			gr->plainfile = NULL;
		}
		else if (!strcmp(query_filename, "-")) {
			gr->plainfile = stdin;
			gr->gzfile = NULL;
			gr->filename = query_filename;
			gr->gz = 0;
		}
		else {
			gr->filename = query_filename;
			gr->plainfile = fopen(gr->filename, "ro");
			gr->gzfile = NULL;
			gr->gz = 0;
		}
	}
	else {
		/*
		 * query is not a valid file name nor standard input, so we return a
		 * NULL GTF_READER
		 */
		free(gr);
		gr = NULL;
	}

	free(tmp);
	return gr;
}


char *readline(TEXTFILE_READER *gr) {
	char *buffer = (char *)calloc(10000, sizeof(char));
	char *ret = NULL, *ret_tmp, *ret_trim;;

	if (gr->gz)
		ret_tmp = gzgets(gr->gzfile, buffer, 10000);
	else
		ret_tmp = fgets(buffer, 10000, gr->plainfile);
	if (ret_tmp != NULL) {
		*(buffer + strlen(buffer) - 1) = 0;
		ret_trim = trim_ip(buffer);
		ret = strdup(ret_trim);
	}
	free(buffer);
	return ret;
}

char *get_blast_header(TEXTFILE_READER *gr, BLAST_HEADER *bh) {
	char *buffer = NULL;
	char **field, *p, *i;

	buffer = readline(gr);
	if (buffer != NULL) {
		bh->program_name = buffer;
		while ((buffer = readline(gr)) != NULL)
			if (strncmp("Query=", buffer, 6)) {
				if (!strncmp("Database:", buffer, 9)) {
					bh->database_name = strdup(strchr(buffer, ' ') + 1);
					free(buffer);
					buffer = readline(gr);
					split_ip(&field, buffer, " ");
					bh->database_nb_sequences = atoi(field[0]);
					i = p = field[2];
					while (*p != 0) {
						if (*p != ',') {
							*i = *p;
							i++;
						}
						p++;
					}
					*i = 0;
					bh->database_length = atoi(field[2]);
					free(field);
					free(buffer);
				}
				else {
					free(buffer);
					buffer = NULL;
				}
			}
			else
				break;
	}
	return buffer;
}

char *get_next_blast_hsp(TEXTFILE_READER *gr, BLAST_HSP *hsp, char *buffer) {
	char *tmp = NULL;
	char *i;
	char **field;
	int deja, nf, k;

	if (buffer != NULL) {
		tmp = strdup(buffer);
		free(buffer);
	}
	else
		tmp = readline(gr);
	deja = 0;

	while (tmp != NULL) {
		if (deja && (!strncmp("Query=", tmp, 6) || (*tmp == '>') || !strncmp("Score", tmp, 5))) {
			break;
		}
		else if (!strncmp("Query=", tmp, 6)) {
			if (hsp->bq.query_name != NULL) free(hsp->bq.query_name);
			hsp->bq.query_name = strdup(strchr(tmp, ' ') + 1);
			free(tmp);
			tmp = readline(gr);
			*strchr(tmp, ' ') = 0;
			hsp->bq.query_length = atoi(tmp + 1);
		}
		else if (*tmp == '>') {
			if (hsp->bs.subject_name != NULL) free(hsp->bs.subject_name);
			hsp->bs.subject_name = strdup(tmp + 1);
		}
		else if (!strncmp("Length", tmp, 6)) {
			i = strchr(tmp, '=');
			if (i != NULL) hsp->bs.subject_length = atoi(trim_ip(i + 1));
		}
		else if (!strncmp("Score =", tmp, 7)) {
			i = strchr(tmp, '=');
			if (i != NULL) {
				split_ip(&field, tmp, " =,");
				hsp->score = atof(field[1]);
				hsp->expect = atof(field[5]);
				free(field);
			}
			deja = 1;
		}
		else if (!strncmp("Identities", tmp, 10)) {
			i = strchr(tmp, '=');
			if (i != NULL) {
				nf = split_ip(&field, tmp, " =%()");
				if (hsp->identities != NULL) free(hsp->identities);
				hsp->identities = strdup(field[1]);
				hsp->identities_percent = atoi(field[2]);
				if (hsp->gaps != NULL) {
					free(hsp->gaps);
					hsp->gaps = NULL;
				}
				if (nf > 3) {
					hsp->gaps = strdup(field[5]);
					hsp->gap_percent = atoi(field[6]);
				}
				free(field);
				hsp->query_start = -1;
				hsp->query_end = -1;
				hsp->subject_start = -1;
				hsp->subject_end = -1;
			}
		}
		else if (!strncmp("Strand", tmp, 6)) {
			i = strchr(tmp, '=');
			if (i != NULL) {
				split_ip(&field, tmp, " =/");
				hsp->strand_query = strcmp(field[1], "Plus") ? '-' : '+';
				hsp->strand_subject = strcmp(field[2], "Plus") ? '-' : '+';
				free(field);
			}
		}
		else if (!strncmp("Query:", tmp, 6)) {
			split_ip(&field, tmp, " ");
			k = atoi(field[1]);
			if ((hsp->query_start == -1) || (hsp->query_start > k))	hsp->query_start = k;
			if ((hsp->query_end == -1) || (hsp->query_end < k))	hsp->query_end = k;
			k = atoi(field[3]);
			if ((hsp->query_start == -1) || (hsp->query_start > k))	hsp->query_start = k;
			if ((hsp->query_end == -1) || (hsp->query_end < k))	hsp->query_end = k;
			free(field);
		}
		else if (!strncmp("Sbjct:", tmp, 6)) {
			split_ip(&field, tmp, " ");
			k = atoi(field[1]);
			if ((hsp->subject_start == -1) || (hsp->subject_start > k))	hsp->subject_start = k;
			if ((hsp->subject_end == -1) || (hsp->subject_end < k))	hsp->subject_end = k;
			k = atoi(field[3]);
			if ((hsp->subject_start == -1) || (hsp->subject_start > k))	hsp->subject_start = k;
			if ((hsp->subject_end == -1) || (hsp->subject_end < k))	hsp->subject_end = k;
			free(field);
		}
		free(tmp);
		tmp = readline(gr);
	}
	return tmp;
}

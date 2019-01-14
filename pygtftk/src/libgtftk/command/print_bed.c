/*
 * write_bed.c
 *
 *  Created on: December 7, 2018
 *      Author: puthier (based on Fafa code...)
 *  Objective: print a gtf obj in bed format
 */

#include "libgtftk.h"

extern void print_row_bed(FILE *output, GTF_ROW *r, char delim, int add_chr, char *keys, char *sep, char *more_info);

__attribute__ ((visibility ("default")))
void *print_bed(GTF_DATA *gtf_data, char *output, int add_chr, char *keys, char *sep, char *more_info) {
	int i;
	FILE *out = stdout;

	if (gtf_data != NULL) {
		if (*output != '-') out = fopen(output, "w");
		if (out == NULL) out = stdout;
		for (i = 0; i < gtf_data->size; i++){
		    print_row_bed(out, gtf_data->data[i], '\t', add_chr, keys, sep, more_info);
		}
		if (out != stdout) {
			fflush(out);
			fclose(out);
		}
	}
	return 0;
}

/*
 * add_prefix.c
 *
 *  Created on: Nov 14, 2017
 *      Author: dputhier (inspired from fafa code...)
 *
 *
 */

#include "libgtftk.h"

/*
 * external functions declaration
 */
extern GTF_DATA *clone_gtf_data(GTF_DATA *gtf_data);
//extern int update_attribute_table(GTF_ROW * row);

/*
 * global variables declaration
 */
extern COLUMN **column;
extern int nb_column;

__attribute__ ((visibility ("default")))
GTF_DATA *add_prefix(GTF_DATA *gtf_data, char *features, char *key, char *txt, int suffix) {
	int i, k, ok;
	int target_field = -1;
	char *str_concat;

	if (strcmp(key, "chrom") == 0) key = "seqid";

	/*
	 * reserve memory for the GTF_DATA structure to return
	 */
	GTF_DATA *ret = clone_gtf_data(gtf_data);

	GTF_ROW *row;
	ATTRIBUTE *pattr;

	/*
	 * Check whether we want to add a prefix/suffix to a basic attribute
	 */
	for (i = 0; i < nb_column - 1; i++)
	    if (!strcmp(column[i]->name, key)) {
	    	target_field = i;
	    	break;
	    }

	for (i = 0; i < ret->size; i++) {
		row = ret->data[i];
		ok = (*features == '*');
		if (!ok) ok = (strstr(features, row->field[2]) != NULL);
		if (ok) {
			/*
			 * if target_field contains -1 will will check extended attributes
			 * else we will modify basic attributes.
			 */
			if (target_field >= 0){
				str_concat = (char *)calloc(strlen(txt) + strlen(row->field[target_field]) + 1, sizeof(char));
				if (suffix) {
					strcpy(str_concat, row->field[target_field]);
					strcat(str_concat, txt);
				}
				else {
					strcpy(str_concat, txt);
					strcat(str_concat, row->field[target_field]);
				}
				free(row->field[target_field]);
				row->field[target_field] = str_concat;
			}
			else {
				/*
				pattr = row->attributes.attr[0];
				while (pattr != NULL) {
					if (strstr(key, pattr->key)) {
						str_concat = (char *)calloc(strlen(txt) + strlen(pattr->value) + 1, sizeof(char));
						if (suffix) {
							strcpy(str_concat, pattr->value);
							strcat(str_concat, txt);
						}
						else {
							strcpy(str_concat, txt);
							strcat(str_concat, pattr->value);
						}
						pattr->value = str_concat;
					}
					pattr = pattr->next;
				}
				*/
				for (k = 0; k < row->attributes.nb; k++) {
					pattr = row->attributes.attr + k;
					if (strstr(key, pattr->key)) {
						str_concat = (char *)calloc(strlen(txt) + strlen(pattr->value) + 1, sizeof(char));
						if (suffix) {
							strcpy(str_concat, pattr->value);
							strcat(str_concat, txt);
						}
						else {
							strcpy(str_concat, txt);
							strcat(str_concat, pattr->value);
						}
						free(pattr->value);
						pattr->value = str_concat;
					}
				}
			}
		}
		//update_attribute_table(row);
	}
	return ret;
}


/*
 * del_attributes.c
 *
 *  Created on: Jun 8, 2017
 *      Author: fafa
 */

#include "libgtftk.h"

/*
 * external functions declaration
 */
extern GTF_DATA *clone_gtf_data(GTF_DATA *gtf_data);
extern int update_attribute_table(GTF_ROW * row);

__attribute__ ((visibility ("default")))
GTF_DATA *del_attributes(GTF_DATA *gtf_data, char *features, char *keys) {
	int i, ok;

	/*
	 * reserve memory for the GTF_DATA structure to return
	 */
	GTF_DATA *ret = clone_gtf_data(gtf_data);

	GTF_ROW *row;
	ATTRIBUTE *pattr, *previous_pattr;

	for (i = 0; i < ret->size; i++) {
		row = ret->data[i];
		ok = (*features == '*');
		if (!ok) ok = (strstr(features, row->field[2]) != NULL);
		if (ok) {
			pattr = row->attributes.attr[0];
			previous_pattr = NULL;
			while (pattr != NULL) {
				if (strstr(keys, pattr->key)) {
					free(pattr->key);
					free(pattr->value);
					if (previous_pattr != NULL)
						previous_pattr->next = pattr->next;
					else
						row->attributes.attr[0] = row->attributes.attr[0]->next;
					row->attributes.nb--;
					pattr = pattr->next;
				}
				else {
					previous_pattr = pattr;
					pattr = pattr->next;
				}
			}
			update_attribute_table(row);
		}
	}

	return ret;
}

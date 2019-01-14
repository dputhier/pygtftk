/*
 * del_attributes.c
 *
 *  Created on: Jun 8, 2017
 *      Author: fafa
 */

#include <libgtftk.h>

/*
 * external functions declaration
 */
extern GTF_DATA *clone_gtf_data(GTF_DATA *gtf_data);
//extern int update_attribute_table(GTF_ROW * row);

__attribute__ ((visibility ("default")))
GTF_DATA *del_attributes(GTF_DATA *gtf_data, char *features, char *keys) {
	int i, ok;

	/*
	 * reserve memory for the GTF_DATA structure to return
	 */
	GTF_DATA *ret = clone_gtf_data(gtf_data);

	GTF_ROW *row;
	ATTRIBUTE *pattr, /**previous_pattr,*/ *pattr_end;

	for (i = 0; i < ret->size; i++) {
		row = ret->data[i];
		ok = (*features == '*');
		if (!ok) ok = (strstr(features, row->field[2]) != NULL);
		if (ok) {
			pattr = row->attributes.attr;
			pattr_end = pattr + row->attributes.nb;
			while (pattr != pattr_end) {
				if (strstr(keys, pattr->key)) {
					if (pattr + 1 != pattr_end)
						bcopy(pattr + 1, pattr, (pattr_end - pattr - 1) * sizeof(ATTRIBUTE));
					row->attributes.nb--;
					pattr_end--;
				}
				else
					pattr++;
			}
			row->attributes.attr = (ATTRIBUTE *)realloc(row->attributes.attr, row->attributes.nb * sizeof(ATTRIBUTE));
		}
	}

	return ret;
}

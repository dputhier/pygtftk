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
extern int update_attribute_table(GTF_ROW * row);
extern void *bookmem(int nb, int size, char *file, const char *func, int line);

/*
 * global variables declaration
 */
extern COLUMN **column;
extern int nb_column;

__attribute__ ((visibility ("default")))
GTF_DATA *add_prefix(GTF_DATA *gtf_data, char *features, char *key, char *txt, int suffix) {
	int i, ok;
	int target_field = -1;

	if(strcmp(key, "chrom") == 0)
		key = "seqid";
	/*
	 * reserve memory for the GTF_DATA structure to return
	 */

	GTF_DATA *ret = clone_gtf_data(gtf_data);

	GTF_ROW *row;

	ATTRIBUTE *pattr;

	/* Check whether we want to add a prefix/suffix to a basic attribute */

//<<<<<<< HEAD
	for (i = 0; i < nb_column - 1; i++) {
	//for(i = 0; i < 8; ++i){
	    if(!strcmp(/*basic_attr[i]*/ column[i]->name, key))
//=======
//	for (i = 0; i < nb_column; i++) {
//
//	    if(!strcmp(column[i]->name, key))
//>>>>>>> c2d0f09934535ced750c90f997c24b68c0686fbc
	    {
	    	target_field = i;
	    	break;
	    }
	}

	for (i = 0; i < ret->size; i++) {
		row = ret->data[i];
		ok = (*features == '*');
		if (!ok) ok = (strstr(features, row->field[2]) != NULL);
		if (ok) {

			/* if target_field contains -1 will will check extended attributes
			 * else we will modify basic attributes.
			 * */
			if(target_field >= 0){
				//char *str_concat = malloc(strlen(txt) + strlen(row->field[target_field]) + 1);
				char *str_concat = (char *)bookmem(strlen(txt) + strlen(row->field[target_field]) + 1, 1, __FILE__, __func__, __LINE__);
				if(suffix){
					strcpy(str_concat, row->field[target_field]);
					strcat(str_concat, txt);
				}else{
					strcpy(str_concat, txt);
					strcat(str_concat, row->field[target_field]);
				}
				row->field[target_field] = str_concat;
			}else{

					pattr = row->attributes.attr[0];
					while (pattr != NULL) {
						if (strstr(key, pattr->key)) {
							//char *str_concat = malloc(strlen(txt) + strlen(pattr->value) + 1);
							char *str_concat = (char *)bookmem(strlen(txt) + strlen(pattr->value) + 1, 1, __FILE__, __func__, __LINE__);
							if(suffix){
								strcpy(str_concat, pattr->value);
								strcat(str_concat, txt);
							}else{
								strcpy(str_concat, txt);
								strcat(str_concat, pattr->value);
							}
							pattr->value = str_concat;
						}
						pattr = pattr->next;
					}
			}
		}

		update_attribute_table(row);
	}
	return ret;

}


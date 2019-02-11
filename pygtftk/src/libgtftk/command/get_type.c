/*
 * get_type.c
 *
 *  Created on: Dec 21, 2018
 *      Author: fafa
 */

#include "libgtftk.h"

extern TTEXT *vret;
extern int nb_column;
extern COLUMN **column;
extern INDEX_ID *index_gtf(GTF_DATA *gtf_data, char *key);

int is_int(char *value) {
	char *p = value;

	if (*p == 0) return 0;
	if (*p == '-' || *p == '+') p++;
	if (*p == 0) return 0;

	while (*p) {
		if (!isdigit(*p)) return 0;
		p++;
	}
	return 1;
}

int is_float(char *value) {
	char *p;
	strtof(value, &p);
	return *p == '\0';
}

int is_boolean(char *value) {
	return (*value == '0' || *value == '1') && *(value + 1) == 0;
}

int is_char(char *value) {
	return isprint((int)*value) && *(value + 1) == '\0';
}

int find_type(char *value) {
	int type;

	if (*value == 0)
		type = IRREGULAR;
	else if ((*value == '.') && (*(value + 1) == 0))
		type = UNAVAILABLE;
	else if ((*value == '?') && (*(value + 1) == 0))
		type = UNKNOWN;
	else if (is_boolean(value))
		type = BOOLEAN;
	else if (is_int(value))
		type = INTEGER;
	else if (is_float(value))
		type = FLOAT;
	else if (is_char(value))
		type = CHAR;
	else
		type = STRING;
	return type;
}

static void action(const void *nodep, const VISIT which, const int depth) {
	ROW_LIST *datap;
	char tmp[100];

	switch (which) {
		case preorder:
			break;
		case postorder:
		case leaf:
			datap = *((ROW_LIST **)nodep);
			vret->data = (char ***)realloc(vret->data, (vret->size + 1) * sizeof(char **));
			vret->data[vret->size] = (char **)calloc(3, sizeof(char *));
			sprintf(tmp, "%d", datap->nb_row);
			vret->data[vret->size][0] = strdup(tmp);
			vret->data[vret->size][1] = strdup(datap->token);
			sprintf(tmp, "%d", find_type(datap->token));
			vret->data[vret->size][2] = strdup(tmp);
			vret->size++;
			break;
		case endorder:
			break;
	}
}

__attribute__ ((visibility ("default")))
int get_type(GTF_DATA *gtf_data, char *key, int ignore_undef) {
	int i, type = UNAVAILABLE, type2;
	INDEX_ID *index_id;
	vret = (TTEXT *)calloc(1, sizeof(TTEXT));

	for (i = 0; i < nb_column - 1; i++)
	    if (!strcmp(column[i]->name, key)) {
	    	index_id = index_gtf(gtf_data, key);
	    	twalk(column[index_id->column]->index[index_id->index_rank]->data, action);
	    	break;
	    }
	if (i == (nb_column - 1)) {
		index_id = index_gtf(gtf_data, key);
		twalk(column[index_id->column]->index[index_id->index_rank]->data, action);
	}

	i = 0;
	type = UNSET;
	while (i < vret->size) {
		if (type == UNSET) {
			type = atoi(vret->data[i][2]);
			if (type == IRREGULAR) break;
		}
		else {
			type2 = atoi(vret->data[i][2]);
			if (type2 == IRREGULAR) {
				type = IRREGULAR;
				break;
			}
			else if (ignore_undef) {
				if (type2 > 0) {
					if (type > 0 && type != type2) {
						type = UNCONSISTENT;
						break;
					}
					else
						type = type2;
				}
			}
			else {
				if (type <= 0 || type2 <= 0) {
					type = UNCONSISTENT;
					break;
				}
				else if (type != type2) {
					type = UNCONSISTENT;
					break;
				}
				else
					type = type2;
			}
		}
		i++;
	}
	return type;
}

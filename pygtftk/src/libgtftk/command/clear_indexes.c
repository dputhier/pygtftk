/*
 * clear_indexes.c
 *
 *  Created on: Jan 30, 2018
 *      Author: fafa
 */
#include "libgtftk.h"

/*
 * global variables in libgtftk.c
 */
extern COLUMN **column;

__attribute__ ((visibility ("default")))

void clear_indexes() {
	int i;

	for (i = 0; i < 9; i++) {
		fprintf(stderr, "%s : %d\n", column[i]->name, column[i]->nb_index);
	}
}

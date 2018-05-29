/*
 * free_raw_data.c
 *
 *  Created on: Apr 4, 2017
 *      Author: fafa
 */

#include "libgtftk.h"

extern void freemem(void *ptr, char *file, const char *func, int line);

__attribute__ ((visibility ("default")))
int free_raw_data(RAW_DATA *raw_data) {
	int i, j;

	//for (i = 0; i < raw_data->nb_columns; i++) free(raw_data->column_name[i]);
	freemem(raw_data->column_name, __FILE__, __func__, __LINE__);
	for (i = 0; i < raw_data->nb_rows; i++) {
		for (j = 0; j < raw_data->nb_columns; j++)
			freemem(raw_data->data[i][j], __FILE__, __func__, __LINE__);
		freemem(raw_data->data[i], __FILE__, __func__, __LINE__);
	}
	freemem(raw_data->data, __FILE__, __func__, __LINE__);
	freemem(raw_data, __FILE__, __func__, __LINE__);
	return 0;
}

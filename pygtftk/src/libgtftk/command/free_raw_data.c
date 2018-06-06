/*
 * free_raw_data.c
 *
 *  Created on: Apr 4, 2017
 *      Author: fafa
 */

#include "libgtftk.h"

__attribute__ ((visibility ("default")))
int free_raw_data(RAW_DATA *raw_data) {
	int i, j;

	for (i = 0; i < raw_data->nb_columns; i++) free(raw_data->column_name[i]);
	free(raw_data->column_name);
	for (i = 0; i < raw_data->nb_rows; i++) {
		for (j = 0; j < raw_data->nb_columns; j++)
			free(raw_data->data[i][j]);
		free(raw_data->data[i]);
	}
	free(raw_data->data);
	free(raw_data);
	return 0;
}

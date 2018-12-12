/*
 * column.c
 *
 *  Created on: Jan 9, 2017
 *      Author: fafa
 *
 *  This source file contains all that is related to column management :
 *  	- print_row to print a whole row with a given delimiter
 *		- index_row to add a row in an index corresponding to the given value
 *		- make_columns/make_column to build the GTF column model
 */

#include "libgtftk.h"

/*
 * external functions declaration
 */
extern int compare_row_list(const void *p1, const void *p2);
extern int split_ip(char ***tab, char *s, char *delim);

/*
 * global variables declaration
 */
extern int nb_column;
extern COLUMN **column;

/*
 * prints a column value (a feature, a start ...) with delimiter.
 * Default values are printed as "."
 *
 * Parameters:
 * 		token:	the character string value to print
 * 		output:	where to print (a file or stdout)
 * 		col:	the corresponding column
 * 		delim:	the delimiter character
 */
void print_string(char *token, FILE *output, COLUMN *col, char delim) {
	if (col->default_value != NULL)
		if (!strcmp(token, col->default_value))
			delim != 0 ? fprintf(output, ".%c", delim) : fprintf(output, ".");
		else
			delim != 0 ? fprintf(output, "%s%c", token, delim) : fprintf(output, "%s", token);
	else
		delim != 0 ? fprintf(output, "%s%c", token, delim) : fprintf(output, "%s", token);
}

/*
 * prints an attribute list from a GTF row.
 *
 * Parameters:
 * 		row:	the GTF row from which to print attributes
 * 		output:	where to print (a file or stdout)
 */
void print_attributes(GTF_ROW *row, FILE *output) {
	int k;
	if (row->attributes.nb != -1) {
		fprintf(output, "%s \"%s\";", row->attributes.attr[0]->key, row->attributes.attr[0]->value);
		for (k = 1; k < row->attributes.nb; k++)
			fprintf(output, " %s \"%s\";", row->attributes.attr[k]->key, row->attributes.attr[k]->value);
	}
}

/*
 * prints a whole GTF row in output with the given character delimiter
 *
 * Parameters:
 * 		output:		where to print
 * 		r:			the row to print
 * 		delim:		the delimiter character
 * 		add_chr:	boolean; if true(1), add "chr" at the begining of the row
 */
void print_row(FILE *output, GTF_ROW *r, char delim, int add_chr) {
	int i;

	if (add_chr) fprintf(output, "chr");
	for (i = 0; i < 8; i++) print_string(r->field[i], output, column[i], delim);
	print_attributes(r, output);
	fprintf(output, "\n");
}

/*
 * prints a list of comma (or other sep) separated attribute values from a GTF row.
 *
 *
 * Parameters:
 * 		row:	the GTF row from which to print attributes
 * 		output:	where to print (a file or stdout)
 *      keys: a comma separated list of keys to be selected.
 *      sep:    the separator to use (e.g '|', ','...).
 */
void print_attr_with_sep(GTF_ROW *row, FILE *output, char delim, char *keys, char *sep, char *more_info) {
	int n, i, j;
	int found = 0;
	char **token;
	char *keys_concat = (char *)calloc(5000, sizeof(char));


    n = split_ip(&token, strdup(keys), ",");

    for(i=0; i<n; i++)
    {
        found = 0;

        if(strcmp("seqid", token[i]) == 0)
        {
                    strcat(keys_concat, row->field[0]);
                    if (i < (n-1)) strcat(keys_concat, sep);
                    found = 1;
        }
        else if(strcmp("chrom", token[i]) == 0)
        {
                    strcat(keys_concat, row->field[0]);
                    if (i < (n-1)) strcat(keys_concat, sep);
                    found = 1;
        }
        else if(strcmp("source", token[i]) == 0)
        {
                    strcat(keys_concat, row->field[1]);
                    if (i < (n-1)) strcat(keys_concat, sep);
                    found = 1;
        }
        else if(strcmp("feature", token[i]) == 0)
        {
                    strcat(keys_concat, row->field[2]);
                    if (i < (n-1)) strcat(keys_concat, sep);
                    found = 1;
        }
        else if(strcmp("start", token[i]) == 0)
        {
                    strcat(keys_concat, row->field[3]);
                    if (i < (n-1)) strcat(keys_concat, sep);
                    found = 1;
        }
        else if(strcmp("end", token[i]) == 0)
        {
                    strcat(keys_concat, row->field[4]);
                    if (i < (n-1)) strcat(keys_concat, sep);
                    found = 1;
        }
        else if(strcmp("score", token[i]) == 0)
        {
                    strcat(keys_concat, row->field[5]);
                    if (i < (n-1)) strcat(keys_concat, sep);
                    found = 1;
        }
        else if(strcmp("strand", token[i]) == 0)
        {
                    strcat(keys_concat, row->field[6]);
                    if (i < (n-1)) strcat(keys_concat, sep);
                    found = 1;
        }
        else if(strcmp("phase", token[i]) == 0)
        {
                    strcat(keys_concat, row->field[7]);
                    if (i < (n-1)) strcat(keys_concat, sep);
                    found = 1;
        }

        if(found == 0)
        {
            if (row->attributes.nb != -1)
            {
                for(j=0; j< row->attributes.nb; j++)
                {

                    if(strcmp(row->attributes.attr[j]->key, token[i]) == 0)
                    {
                            strcat(keys_concat, row->attributes.attr[j]->value);
                            if (i < (n-1)) strcat(keys_concat, sep);
                            found = 1;

                    }
                }
            }
        }

        if(found == 0){
            if(found == 0){
                strcat(keys_concat, "?");
                if (i < (n-1)) strcat(keys_concat, sep);
            }
        }
    }

    //keys_concat = (char *)realloc(keys_concat, (strlen(keys_concat) + 1) * sizeof(char));


    if(strlen(more_info)){
        strcat(keys_concat, sep);
        strcat(keys_concat, more_info);
    }

    if(keys_concat){
       fprintf(output, "%s", keys_concat);
    }
    else
    {
        fprintf(output, "");
    }

    fprintf(output, "%c", delim);

    free(keys_concat);
}

/*
 * prints a GTF row in bed format to output with the given character delimiter
 *
 * Parameters:
 * 		output:		where to print
 * 		r:			the row to print
 * 		delim:		the delimiter character
 * 		add_chr:	boolean; if true(1), add "chr" at the begining of the row
 *      keys:      a comma separated list of keys whose values will appear in the name
 *                  (4th) colum of the bed.
 *      sep:        The separator to use for names (4th column).
 */
void print_row_bed(FILE *output, GTF_ROW *r,  char delim, int add_chr, char *keys, char *sep, char *more_info){


	if (add_chr) fprintf(output, "chr");
	// print chr, start, end
	print_string(r->field[0], output, column[0], delim);
	fprintf(output,"%d%c", atoi(r->field[3]) - 1, delim);
	print_string(r->field[4], output, column[4], delim);
	// print requested columns with user-defined delim
	print_attr_with_sep(r, output, delim, keys, sep, more_info);
	// print score strand
	print_string(r->field[5], output, column[5], delim);
	fprintf(output,"%s", r->field[6]);
	fprintf(output, "\n");
}


/*
 * Creates an empty index and add/link it in the index list of the column
 * in index_id.
 *
 * Parameters:
 * 		index_id:	contains the column to index
 * 		key:		the column name or attribute name to index
 */
void make_index(INDEX_ID *index_id, char *key) {
	column[index_id->column]->nb_index++;
	column[index_id->column]->index = (INDEX **)realloc(column[index_id->column]->index, column[index_id->column]->nb_index * sizeof(INDEX *));
	INDEX *index = (INDEX *)calloc(1, sizeof(INDEX));
	column[index_id->column]->index[column[index_id->column]->nb_index - 1] = index;
	index->data = NULL;
	index->gtf_data = NULL;
	index->key = strdup(key);
	index_id->index_rank = column[index_id->column]->nb_index - 1;
	if (column[index_id->column]->nb_index > 1)
		column[index_id->column]->index[column[index_id->column]->nb_index - 2]->next = index;
}

/*
 * This function indexes a row.
 *
 * Parameters:
 * 		row_nb:		the row number to add
 * 		value:		the value with which to row is associated
 * 		index:		the index of ROW_LIST elements
 */
void index_row(int row_nb, char *value, INDEX *index) {
	ROW_LIST *test_row_list, *row_list, *find_row_list;

	if (index != NULL) {
		/*
		 * build a ROW_LIST to check if value is already indexed
		 */
		test_row_list = calloc(1, sizeof(ROW_LIST));

		test_row_list->token = value;
		find_row_list = tfind(test_row_list, &(index->data), compare_row_list);
		if (find_row_list == NULL) {
			/*
			 * value is not in the index so we reserve a ROW_LIST, initialize it
			 * with one row (row_nb) and put it in the index (tsearch)
			 */

			row_list = (ROW_LIST *)calloc(1, sizeof(ROW_LIST));
			row_list->token = value;
			row_list->nb_row = 1;
			row_list->row = (int *)calloc(1, sizeof(int));
			row_list->row[0] = row_nb;
			tsearch(row_list, &(index->data), compare_row_list);
		}
		else {
			/*
			 * value is already in the index so we just have to add row_nb into
			 * the ROW_LIST element found
			 */
			row_list = *((ROW_LIST **)find_row_list);
			row_list->nb_row++;
			row_list->row = (int *)realloc(row_list->row, row_list->nb_row * sizeof(int));
			row_list->row[row_list->nb_row - 1] = row_nb;
		}
		free(test_row_list);
	}
}

/*
 * Reserve memory for a COLUMN structure and fill it whith parameters.
 * This function doesn't create any index in columns.
 *
 * Parameters:
 * 		i:		the number of the column (0 to 8)
 * 		dv:		the default value
 * 		name:	the column name
 *
 * Returns:		the initialized COLUMN structure
 */
COLUMN *make_column(int i, void *dv, char *name) {
	COLUMN *column = (COLUMN *)calloc(1, sizeof(COLUMN));
	column->num = i;
	column->name = strdup(name);
	column->index = NULL;
	column->nb_index = 0;
	if (dv != NULL) column->default_value = strdup((char *)dv);

	return column;
}

/*
 * Creates and initialize a column model suitable for GTF data.
 */
void make_columns(void) {
	nb_column = 9;
	if (column == NULL) {
		column = (COLUMN **)calloc(nb_column, sizeof(COLUMN *));
		column[0] = make_column(0, ".", "seqid");
		column[1] = make_column(1, ".", "source");
		column[2] = make_column(2, ".", "feature");
		column[3] = make_column(3, ".", "start");
		column[4] = make_column(4, ".", "end");
		column[5] = make_column(5, ".", "score");
		column[6] = make_column(6, ".", "strand");
		column[7] = make_column(7, ".", "phase");
		column[8] = make_column(8, ".", "attributes");
	}
}

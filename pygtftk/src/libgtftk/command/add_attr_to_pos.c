/*
 * add_attr_column.c
 *
 *  Created on: Nov 14, 2017
 *      Author: dputhier.
 *      Purpose: Simply add a new key at specified line
 *      inputfile_name is a two column file with position (line number) as column 1 and values as column 2.
 *      TODO: should check that column 1 is an int otherwise, stop.
 */

#include "libgtftk.h"
#include <ctype.h>

/*
 * external functions declaration
 *
 */
extern GTF_DATA *clone_gtf_data(GTF_DATA *gtf_data);
extern void add_attribute(GTF_ROW *row, char *key, char *value);
extern int split_ip(char ***tab, char *s, char *delim);

/*
 * global variables declaration
 */

__attribute__ ((visibility ("default")))
GTF_DATA *add_attr_to_pos(GTF_DATA *gtf_data, char *inputfile_name, char *new_key) {

    /*
     * Required variables
     */

    int n;
    GTF_DATA *ret = clone_gtf_data(gtf_data);
    GTF_ROW *row;
    FILE *input = fopen(inputfile_name, "ro");
    size_t buffersize = 1000;
    char *buffer = (char *)calloc(buffersize, sizeof(char));
    char **token;


    /*
     * loop on the txt file rows
     */

    while (fgets(buffer, 999, input) != NULL) {

        if(strcmp(buffer,"\n") || strcmp(buffer,"\r\n")){
            // Delete newline
            strtok(buffer, "\r\n");

            /*
             * split the line and check the number of fields (should be 2)
             */
            n = split_ip(&token, buffer, "\t");


            if (n < 2) {
                fprintf(stderr, "ERROR : Need two columns as inputfile (add_attr_to_pos).\n");
                exit(0);
            }

            if (n > 2) {
                fprintf(stderr, "ERROR : need two columns.");
                exit(0);
            }

            if(atoi(token[0]) <= ret->size){
                row = ret->data[atoi(token[0])];
                add_attribute(row, new_key, token[1]);
             }else{
                fprintf(stderr, "ERROR : index out of range (add_attr_to_pos).");
                exit(0);
             }


            }

    }
    free(buffer);
    return ret;
}

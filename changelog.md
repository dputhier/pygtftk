## Development

Enhancements:


* Added a specific object for Exceptions (GTFtkError). This object knows if a message should be printed when gtftk is used or not in interactive mode.
* Included call to add_attributes (libgtftk) into join_attr, convergent, divergent, closest_gene, nb_exons and overlapping commands.
* Added new tests to join_attr command.
* convert_to_ensembl is now propose as a methode for GTF classe. The method use the C version (libgtftk). New tests have been added to convert_ensembl command.
* Introduce -a to gtftk in order to download new plugins. Still experimental.
* The -l argument to gtftk now allows to list available commands.
* A new command wise argument (-C) allow to add chr to the chromosome names when write method is called. 
* All plugins now require a test !
* New tests have been added to heatmap and profile command.
* Added the select_by_numeric_value command.
* Added the apropos command to search through command descriptions.
* Added the select_by_regexp command.
* Added the select_by_position command.
* Added the get_attr_value_list command.
* Added the col_from_tab command.
* Added the control_list command.
* Added the select_by_intron_size command.
* Added the shift command.
* The cloud dependency has been changed by a cloudpickle dependency.
* Added seqid_list.
* Added a --explicit argument to get_tx_seq.
* GTF is now the default output of nb_exons, nb_transcripts, closest_genes.
* Added pyparsing and biopython as dependencies.
* Added new tests.
* The commands add_prefix, merge_attr now rely on C calls through libgtftk.

Bug Fixes:

* Fixed errors and warning in documentation. 


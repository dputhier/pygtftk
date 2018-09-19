# libgtftk README section

Libgtftk is a core library designed to parse and index GTF files. This library provides basic functions to loop through a GTF file efficiently, and to perform basic operations including selecting/subsetting by keys or features, editing/updating key values and more complexe operations at the gene or transcript level such as fetching feature sequences from FASTA files.. This library is intended to be called by external clients written in various languages including C, JAVA, [Python](https://github.com/dputhier/pygtftk), R...  

## System requirements

The libgtftk has been tested on UNIX-like systems including:

 - OSX (10.11.6)
 - Linux (Ubuntu 12.04, 14.04, 16.04)


## Compilation


Under Linux:

      git clone https://github.com/dputhier/libgtftk.git
      cd libgtftk
      gcc -o libgtftk.so -fPIC *.c command/*.c -I. -lz -shared

Under OSX:

      git clone https://github.com/dputhier/libgtftk.git
      cd libgtftk
      gcc -o libgtftk.so -fPIC *.c command/*.c -I. -lz -dynamiclib


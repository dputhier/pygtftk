.. gtftk documentation master file, created by
   sphinx-quickstart on Fri Jan  2 11:18:01 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Gtftk documentation page
========================

The **Python GTF toolkit (pygtftk) package** is intented to ease handling of GTF (Gene Transfer Format) files. The pygtftk package is compatible with Python 2.7 and relies on **libgtftk**, a library of functions **written in C**.

The package comes with a set of **UNIX commands** that can be accessed through the **gtftk  program**. The gtftk program proposes several atomic tools than can be piped to filter, convert, or extract data from GTF files (...). The gtftk set of Unix commands, can be easily extended using a basic plugin architecture. All these aspects are covered in the help section.

While the gtftk Unix program comes with hundreds of unitary and functionnal tests, it is still upon  active development and may thus suffer from bugs that remain to be discovered. Feel free to post any problem or required enhancement in the issue section of the github repository.

.. toctree::
   :maxdepth: 2
   
   installation
   presentation
   api
   bwig_coverage
   developers


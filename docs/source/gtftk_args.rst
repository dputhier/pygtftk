Help on gtftk Unix commands
============================


Main parser arguments of gtftk
-------------------------------


Getting help with -h
~~~~~~~~~~~~~~~~~~~~~

The -h argument can be used to get a synopsis for implemented commands.

.. command-output:: gtftk -h
	:shell:


------------------------------------------------------------------------------------------------------------------


Activating Bash completion
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The code provided below can be useful to activate bash completion.

.. code-block:: bash

  # Use the -b argument of gtftk
  # This will produce a script that you
  # should store in your .bashrc
  gtftk -b

Or alternatively

.. code-block:: bash

  echo "" >> ~/.bashrc
  gtftk -b >> ~/.bashrc



------------------------------------------------------------------------------------------------------------------


Getting the list of funtional tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One can access to the list of implemented  tests through the -p/--plugin-tests arguments. These tests may be run using `bats <https://github.com/sstephenson/bats>`_ (Bash Automated Testing System).


.. code-block:: bash

  # gtftk --plugin-tests



------------------------------------------------------------------------------------------------------------------


Command-wide arguments
--------------------------

**Description:** The following arguments are available in almost all gtftk commands :

- -h, --help : Argument list and details.
- -i, --inputfile: The input file (may be <stdin>).
- -o, --outputfile: The output file (may be <stdout>).
- -D, --no-date: Do not add date to output file names.
- -C, --add-chr: Add 'chr' to chromosome names before printing output.
- -V, --verbosity: Increases output verbosity (can take value from 0 to 4).
- -K --tmp-dir: Keep all temporary files into this folder.
- -L, --logger-file: Stores the values of all command line arguments into a file.


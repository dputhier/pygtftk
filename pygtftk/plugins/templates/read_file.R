#!/usr/bin/env RScript
version = 0.1
desc <- paste("A plugin in R that print something.",
			  "Developer: D. Puthier.")


suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description=desc,
						conflict_handler='resolve')

parser$add_argument('-f',
					'--infile',  
					metavar='\b',  
					default="",
					help='a character')

parser$add_argument('-v',
					'--verbosity',  
					action="store_true",  
					help='Increase output verbosity. [default 0]')
			
parser$add_argument('-i',
					'--aninteger',  
					type='integer',
					metavar='\b',
					default='1',
					help='An integer')

parser$add_argument('-o', '--output',  
					type="integer",
					metavar='\b',
					help='Increase output verbosity. [default 0]',
					required=TRUE)
			

# This will be true unless
# a single line with "options(run.main=FALSE)"
# is added before sourcing this code.

if (getOption('run.main', default=TRUE)) {
	parser$parse_args()
	# Do something
	cat("Something that remains to be defined...\n")
}



# ROSE : RANK ORDERING OF SUPER-ENHANCERS

CLONED using SOURCETREE from: https://bitbucket.org/young_computation/rose/src/master/

### Changes from Source.
1. USAGE

	```bash
	PATHTO=/path/to/ROSE
	PYTHONPATH=$PATHTO/lib
	export PYTHONPATH
	export PATH=$PATH:$PATHTO/bin

	ROSE_main.py [options] -g [GENOME] -i [INPUT_REGION_GFF] -r [RANKBY_BAM_FILE] -o [OUTPUT_FOLDER] [OPTIONAL_FLAGS]
	```

1. Update: 

	* ROSE is executable independent of software directory location.
	* ROSE is compatible with python3
	* ROSE as a docker image ghcr.io/stjude/rose:latest

1. REQUIREMENTS:

	1. All files :
	All input files much be in one directory.

	1. Annotation file :
	Annotation file should be in UCSC table track format (https://genome.ucsc.edu/cgi-bin/hgTables).
	Annotation file should be saved as [GENOME]_refseq.ucsc (example: hg19_refseq.ucsc).
	Annotation file should be in annotation/ folder in the input files directory.

	1. BAM files (of sequencing reads for factor of interest and control) :
	Files must have chromosome IDs starting with "chr"
	Files must be sorted and indexed using SAMtools in order for bamToGFF.py to work. (http://samtools.sourceforge.net/samtools.shtml)

	1. Peak file of constituent enhancers :
	File must be in GFF format with the following columns:

			column 1: chromosome (chr#)
			column 2: unique ID for each constituent enhancer region
			column 4: start of constituent
			column 5: end of constituent
			column 7: strand (+,-,.)
			column 9: unique ID for each constituent enhancer region
			
		NOTE: if value for column 2 and 9 differ, value in column 2 will be used

1. DIRECTORY structure
	```
	├── LICENSE.txt
	│
	├── README.md
	│
	├── bin
	│   ├── ROSE_bamToGFF.py : calculates density of .bam reads in .gff regions
	│   ├── ROSE_callSuper.R : ranks regions by their densities, creates cutoff
	│   ├── ROSE_geneMapper.py : assigns stitched enhancers to genes
	│   └── ROSE_main.py : main program
	└── lib
	    └── ROSE_utils.py : utilities method

	Total: 2 directories, 8 files
	```
1. DEPENDENCIES

	* samtools
	* R version > 3.4
	* bedtools > 2
	* python3

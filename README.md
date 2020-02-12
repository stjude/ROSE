===== ROSE: RANK ORDERING OF SUPER-ENHANCERS =====

LICENSE: Information in LICENSE.txt. 

CLONED using SOURCETREE from: https://bitbucket.org/young_computation/rose/src/master/


#### === Changelog
1. USAGE

	- Option 1: To run the program locally and independent of software location by calling ROSE-local.sh
		NB: First open ROSE-local.sh and modify PATHTO with the PATH ROSE is installed in.
		> ROSE-local.sh ["GTF file"] ["BAM file"] ["OutputDir"] ["feature type"] ["species"] ["bed fileA"] ["bed fileB"]

	- Option 2: Add ROSE to user executable $PATH
		```bash
		PATHTO=/path/to/ROSE
		PYTHONPATH=$PATHTO/lib
		export PYTHONPATH
		export PATH=$PATH:$PATHTO/bin
		```

1. Update: 

	* ROSE is executable independent of software directory location.
	* ROSE has a wrapper script "ROSE-local.sh" to successfully execute all steps of the package, else add ROSE to user executable $PATH
	* ROSE is compatible with python3

1. DIRECTORY structure
	```
	├── LICENSE.txt
	│
	├── README.md
	│
	├── ROSE-local.sh : bash wrapper
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

=================================================================

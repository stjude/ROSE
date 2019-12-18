===== ROSE: RANK ORDERING OF SUPER-ENHANCERS =====

LICENSE: Information in LICENSE.txt. 

MODIFIED BY :

		Modupeore Adetunji
		Abraham's Lab. Computational Biology Department.
		St. Jude Children's Research Hospital.


#### === Changelog 
1) DIRECTORY structure

    ├── LICENSE.txt

    ├── README.md

    ├── ROSE-call.sh    : bash wrapper

    ├── lib

        └── ROSE_utils.py   : utilities method
    │   

    └── src
    
        ├── ROSE_bamToGFF.py    : calculates density of .bam reads in .gff regions
    
        ├── ROSE_callSuper.R    : ranks regions by their densities, creates cutoff
    
        ├── ROSE_geneMapper.py  : assigns stitched enhancers to genes
    
        └── ROSE_main.py    : main program

    2 directories, 8 files

2) DEPENDENCIES

* samtools
* R version > 3.4
* python2


3) USAGE

Program is run by calling ROSE-call.sh
NB: First modify PATHTO with your own directory. 
> ROSE-call.sh ["GTF file"] ["BAM file"] ["OutputDir"] ["feature type"] ["species"] ["bed fileA"] ["bed fileB"]

=================================================================

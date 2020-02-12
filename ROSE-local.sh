#!/bin/bash
#
# Rose Caller to detect both Enhancers and Super-Enhancers
#
# Version 1 11/16/2019

##############################################################
# ##### Please replace PATHTO with your own directory ###### #
##############################################################
PATHTO=/path/to/ROSE
PYTHONPATH=$PATHTO/lib
export PYTHONPATH
export PATH=$PATH:$PATHTO/bin

if [ $# -lt 7 ]; then
  echo ""
  echo 1>&2 Usage: $0 ["GTF file"] ["BAM file"] ["OutputDir"] ["feature type"] ["species"] ["bed fileA"] ["bed fileB"]
  echo ""
  exit 1
fi

#================================================================================
#Parameters for running

# GTF files
GTFFILE=$1

# BAM file
BAMFILE=$2

# Output Directory
OUTPUTDIR=$3
OUTPUTDIR=${OUTPUTDIR:=ROSE_out}

# Feature type
FEATURE=$4
FEATURE=${FEATURE:=gene}

# Species
SPECIES=$5
SPECIES=${SPECIES:=hg_19}

# Bed File A
FILEA=$6

# Bed File B
FILEB=$7

# Transcription Start Size Window
#TSS=
TSS=${TSS:=2000}

# Maximum linking distance for stitching
#STITCH=
STITCH=${STITCH:=12500}


echo "#############################################"
echo "######             ROSE v1             ######"
echo "#############################################"

echo "Input Bed File A: $FILEA"
echo "Input Bed File B: $FILEB"
echo "BAM file: $BAMFILE"
echo "Output directory: $OUTPUTDIR"
echo "Species: $SPECIES"
echo "Feature type: $FEATURE"
#================================================================================
# 
# GENERATING UCSC REFSEQ FILE
#
mkdir -p annotation
echo -e "#bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\tX\tX\tX\t\tX\tname2" > annotation/$SPECIES"_refseq.ucsc"

if [[ $FEATURE == "gene" ]]; then
awk -F'[\t ]' '{
  if($3=="gene")
    print "0\t" $14 "\tchr" $1 "\t" $7 "\t" $4 "\t" $5 "\t" $4 "\t" $5 "\t.\t.\t.\t.\t" $18}' $GTFFILE | sed s/\"//g >> annotation/$SPECIES"_refseq.ucsc"

elif [[ $FEATURE == "transcript" ]]; then
awk -F'[\t ]' '{
  if($3=="transcript")
    print "0\t" $14 "\tchr" $1 "\t" $7 "\t" $4 "\t" $5 "\t" $4 "\t" $5 "\t.\t.\t.\t.\t" $18}' $GTFFILE | sed s/\"//g >> annotation/$SPECIES"_refseq.ucsc"
fi
echo "Annotation file: "$SPECIES"_refseq.ucsc"

#
# GENERATING MERGE BED FILES
#
cat $FILEA $FILEB | sort -k1,1 -k2,2n | mergeBed -i - | awk -F\\t '{print $1 "\t" NR "\t\t" $2 "\t" $3 "\t\t.\t\t" NR}' > unionpeaks.gff
echo "Merge Bed file: unionpeaks.gff"
echo

#
# ROSE CALLER
#
ROSE_main.py -s $STITCH -t $TSS -g $SPECIES -i unionpeaks.gff -r $BAMFILE -o $OUTPUTDIR

echo "Done!"

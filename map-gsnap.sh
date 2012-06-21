#!/bin/bash
FILE=$1
BASENAME=`basename $FILE .fastq`
echo $FILE
GENOME_DIR=data/athaliana-genome/athaliana10/
GENOME=athaliana10
OUTPUT_DIR=data/alignments

echo "aligning $FILE ..."
gsnap -m3 -N1 -A sam -D $GENOME_DIR -d $GENOME $FILE > $OUTPUT_DIR/$BASENAME.sam
echo "completed $FILE"
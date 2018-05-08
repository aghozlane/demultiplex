# Demultiplex
Simple program to demultiplex paired-end and single end reads


## Help
python demultiplex.py -h

usage: demultiplex.py -h

optional arguments:
  -h, --help            show this help message and exit
  -i1 FASTQ_FILE_R1     Fastq R1 file to demultiplex.
  -i2 FASTQ_FILE_R2     Fastq R2 file to demultiplex.
  -i FASTQ_FILE         Single end Fastq file to demultiplex.
  -a ASSIGNATION_TAB_FILE
                        Table assigning index to sample.
  -o OUTPUT_DIR         Output directory (default ./).

Assignation table must be in the tsv format:
barcode sample

#!/bin/bash

query_file=$1
subject_file=$2
output_file=$3

blastn -query $1 -subject $2 -task blastn-short -outfmt 6 | awk '$3==100 && $4==28 {print $0}' > $3 

perfect_matches=`wc -l "$3" | cut -d" " -f1`

echo "Number of perfect matches in $1 is: $perfect_matches"
cat $3

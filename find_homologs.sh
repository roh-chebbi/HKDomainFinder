#!/bin/bash

query_file=$1
subject_file=$2
output_file=$3

tblastn -query $1 -subject $2 -outfmt 6 | awk '{ if ($3 > 30 && $4>0.9*($8-$7)) {print $0} }' 

approx_matches=`wc -l "$3" | cut -d" " -f1`

echo "Number of approx matches in $1 is: $approx_matches"

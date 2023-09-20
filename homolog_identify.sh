#!/bin/bash

query_file=$1
subject_file=$2
output_file=$3

tblastn -query $1 -subject $2 -outfmt "6 std qlen" | awk '{ if ($3 > 30 && $4>0.9*($13)) {print $0} }' > $3 

approx_matches=`wc -l "$3" | cut -d" " -f1`

echo -e "Number of approx matches in $2 is: $approx_matches\n"

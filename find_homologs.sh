#!/bin/bash

query_file=$1
subject_file=$2
output_file=$3

tblastn -query $1 -subject $2 -outfmt 6 | awk '{ if ($3 > 30 && $4>0.9*($8-$7)) {print $0} }' 


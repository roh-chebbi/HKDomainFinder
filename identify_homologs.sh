#!/bin/bash

query_file=$1
subject_file=$2
bed_file=$3
output_file=$4

file_name=$(basename "$2")
org_name="${file_name%.*}"

blast_out=$(mktemp)
genes_out=$(mktemp)

tblastn -query $query_file -subject $subject_file -outfmt "6 std qlen" | awk '$3 > 30 && $4>0.9*$13' \
| cut -f2,9,10 > $blast_out 

while read sid sstart send
do
	awk -F'\t' -v sstart="$sstart" -v send=$send -v sid=$sid \
	'{if ($2<sstart && sstart<$3 && $2<send && send<$3) print $4}' $bed_file >> $genes_out 
done < $blast_out

cat $genes_out | sort | uniq > $output_file
num_genes=$(wc -l < $output_file)

echo -e "Number of unique genes in $org_name is: $num_genes \n"

rm $blast_out
rm $genes_out

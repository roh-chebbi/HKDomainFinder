#!/bin/bash

query_file=$1
subject_file=$2
bed_file=$3
outfile=$4

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

uniq_genes=$(sort $genes_out | uniq)
num_genes=$(echo $uniq_genes | wc -w)

echo -e "Number of genes in $org_name is: $num_genes \n"
echo $org_name $uniq_genes >> $outfile

rm $blast_out
rm $genes_out

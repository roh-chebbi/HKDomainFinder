# HKDomainFinder

This Bash script identifies genes containing homologous histidine kinase (HK) domains in a subject genome. It uses `tblastn` for alignments and checks if these domains are fully contained within gene boundaries from a BED file.

## Requirements
- **BLAST+**: Install via Conda using `conda install -c bioconda blast`.
- **HK Domain faa file**
- **FNA file for each organism**
- **BED file for each organism**

## Usage
```bash
./homolog_identify.sh <query.faa> <subject.fna> <bedfile>
<outfile>  
```

## Example
```bash
for i in {Wolbachia,Vibrio_cholerae_N16961,Escherichia_coli_K12}; do ./homolog_identify.sh data/HK_domain.faa data/"$i".fna data/"$i".bed "$i"_output.txt; done
```

## Outputs
- **Unique Gene Names**: Written to `<output.txt>.`.
- **Console Summary**: Displays the number of unique HK-containing genes.

<img width="1325" alt="image" src="https://github.com/user-attachments/assets/8bf98667-8921-4f7d-bc2c-7eb3b7564a74" />


## How It Works
1. **BLAST Alignment**: `tblastn` aligns query proteins against the subject genome (30% sequence identity & 90% query coverage)
2. **Domain-Gene Mapping**: Checks if BLAST hits are within gene boundaries (gene_start < domain_start < domain_end < gene_end)
3. **Output Generation**: Writes unique gene names to the output file.

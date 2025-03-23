# HKDomainFinder

This tool identifies genes containing homologous histidine kinase (HK) domains in a subject genome. It uses `tblastn` for alignments and checks if these domains are fully contained within gene boundaries from a BED file.

Available in `Bash` and `Python` implementations

## Common Requirements
- **BLAST+**: Install via Conda using `conda install -c bioconda blast`.
- **HK Domain faa file**
- **FNA file for each organism**
- **BED file for each organism**
### Python Specific Requirements
- pandas: pip install pandas or conda install pandas
- subprocess: pip install subprocess or conda install subprocess



## Usage
Bash
```bash
./homolog_identify.sh <query.faa> <subject.fna> <bedfile> <outfile>  
```

Python

```python
python homology_identify.py <blast output file> <bed file> <assembly file> <output file>
```

## Example

Batch Processing in Bash
```bash
for i in {Wolbachia,Vibrio_cholerae_N16961,Escherichia_coli_K12}; do ./homolog_identify.sh data/HK_domain.faa data/"$i".fna data/"$i".bed "$i"_output.txt; done
```
Python
```python
python identify_homologs.py blast_output.txt data/HK_domain.faa data/Vibrio_cholerae_N16961.bed data/Vibrio_cholerae_N16961.fna output.txt
```

## Outputs
- **Unique Gene Names**: Written to `<output.txt>.`.
- **Console Summary**: Displays the number of unique HK-containing genes.
- Python version additionally outputs:

    FASTA File: Sequence records for matched genes


<img width="1325" alt="image" src="https://github.com/user-attachments/assets/8bf98667-8921-4f7d-bc2c-7eb3b7564a74" />


## How It Works
1. **BLAST Alignment**: `tblastn` aligns query proteins against the subject genome (30% sequence identity & 90% query coverage)
2. **Domain-Gene Mapping**: Checks if BLAST hits are within gene boundaries (gene_start < domain_start < domain_end < gene_end)
3. **Output Generation**: Writes unique gene names to the output file.
4. Python Implementation

    BLAST result parsing with pandas

    Biopython-based sequence extraction

    Reverse complement for - strand genes

    FASTA-formatted output


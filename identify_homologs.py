#!/usr/bin/env python

import sys
import pandas as pd
import subprocess

import subprocess

def run_tblastn(query, subject, output_file):
    cmd = [
        "tblastn",
        "-query", query,
        "-subject", subject,
        "-outfmt", "6 std qlen",
        "-out", output_file
    ]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running BLAST: {e}")
        raise
    
    return output_file


def parse_assembly(assembly_file):
    assembly_dict={}
    f = open(assembly_file, "r")

    for line in f:
        if line.startswith(">"):
            key=line.split(">")[1].split(" ")[0]
            assembly_string=[]
        if not line.startswith(">"):
            assembly_string.append(line.rstrip())
        assembly_dict[key]="".join(assembly_string)

    f.close()

    return(assembly_dict)

def filter_blast(blast_file):
    df1 = pd.read_csv(blast_file,sep="\t", header=None)
    df1.columns = ["query acc.ver", "subject acc.ver", "perc identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score", "query length"]

    #filtering based on perc. identity and query length
    df1 = df1[df1["perc identity"]>30]
    df1 = df1[df1["alignment length"]>0.9*df1["query length"]]

    df3=df1[["subject acc.ver","s. start", "s. end"]]
    
    return(df3)

def identify_bed(bed_file,df3):
    df2 = pd.read_csv(bed_file,sep="\t", header=None)
    df5=pd.DataFrame()

    gene_list = []
    for index1,row1 in df3.iterrows():
        sub_start   = row1["s. start"]
        sub_end     = row1["s. end"]
        sub_id      = row1["subject acc.ver"]
        for index2,row2 in df2.iterrows():
            bed_start   = row2[1]
            bed_end     = row2[2]
            gene_name   = row2[3]
            #orientation = row2[5]
            bed_id      = row2[0]
            # check for sequence ID
            if sub_id != bed_id:
                continue
            # check if subject is within bedfile coordinate boundaries
            if (sub_start > bed_start and sub_start < bed_end and sub_end > bed_start and sub_end < bed_end ):
                if not gene_name in gene_list:
                    #check for duplicates
                    gene_list.append(gene_name)
                    df5=pd.concat([df5,row2.to_frame().T])

    df5 = df5.drop(df5.columns[[4]], axis=1)
    df5.columns=["Bed ID","Bed Start","Bed End","Gene Name", "Orientation"]

    return(df5)

def reverse_complement(sequence):
    return(sequence.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])

def extract_sequence(assembly_dict,df,output_file):
    f = open(output_file,"a")
    unique_sequences = set()

    for index,row in df.iterrows():
        bed_start   = int(row["Bed Start"])
        bed_end     = int(row["Bed End"])
        gene_name   = row["Gene Name"]
        orientation = row["Orientation"]
        bed_id      = row["Bed ID"]
        sequence = assembly_dict[bed_id][bed_start-1:bed_end]
        if orientation == "-":
            sequence = reverse_complement(sequence)
        
        unique_sequences.add(sequence)

        f.write(f">{gene_name}\n")
        f.write(f"{sequence}\n")
    
    f.close()
    print(f"Number of Homologs: {len(unique_sequences)} ")

blast_file      = sys.argv[1]
HK_file         = sys.argv[2]
bed_file        = sys.argv[3]
assembly_file   = sys.argv[4]
output_file     = sys.argv[5]

run_tblastn(HK_file, assembly_file,blast_file)
assembly_dict = parse_assembly(assembly_file)
df1 = filter_blast(blast_file)
df2 = identify_bed(bed_file,df1)
extract_sequence(assembly_dict,df2,output_file)
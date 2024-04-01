import json
import csv
import pandas as pd

universal = {
    "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*",
    "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"
    }

def longest_substring_divisible_by_3(s):
    remainder = len(s) % 3 
    length_divisible_by_three = len(s) - remainder 
    return s[:length_divisible_by_three]


def translate(codon):
    
    aa = ''
    codon = codon.upper()
    if len(codon) == 3:
        dictkey=codon.replace("T","U")
        aa+=universal[dictkey]

    return aa

def calculate_cds_coordinates(gene_data):
    cds_coordinates = {}
    position_in_cds = 0  
    cds_exons = gene_data['cds'] if gene_data['strand'] == '+' else reversed(gene_data['cds'])
    for exon in cds_exons:
        coordinates = list(map(int, exon['coordinates']))
        sequence = exon['seq']
        if gene_data['strand'] == '+':
            start_coordinate = coordinates[0]
            for i in range(len(sequence)):
                cds_coordinates[position_in_cds + i + 1] = start_coordinate + i
        else:  # '-' strand
            start_coordinate = coordinates[1]
            for i in range(len(sequence)):
                cds_coordinates[position_in_cds + i + 1] = start_coordinate - i

        position_in_cds += len(sequence)  
    return cds_coordinates

def get_alt(cds, gene_data):
    subs_dict = {}
    cds_coordinates = calculate_cds_coordinates(gene_data)
    for i in range(0, len(cds), 3):
        codon = cds[i:i+3]
        ref_aa = translate(codon)
        if len(codon) == 3:
            for j in range(3):
                ref_nuc = codon[j]
                nuc_coordinate = cds_coordinates[i+j+1]
                for nuc in ["A", "T", "C", "G"]:
                    if ref_nuc != nuc:
                        alt_codon = codon[:j] + nuc + codon[j+1:]
                        alt_aa = translate(alt_codon)
                        aa_pos = (i // 3) + 1
                        key = (i//3, j, ref_nuc, nuc)
                        subs_dict[key] = {
                            'coordinate': nuc_coordinate,
                            'ref_aa': ref_aa,
                            'alt_aa': alt_aa,
                            'aa_pos': aa_pos
                        }
                        
    return subs_dict

file_path = 'all_canonical_transcripts.json'

with open(file_path, 'r') as file:
    data = json.load(file)

all_substitutions_dict = {}

for transcript in data.keys():
    all_substitutions_dict[transcript] = {}
    cds = ""
    exon_coordinates_dict = calculate_cds_coordinates(data[transcript])
    chrm = data[transcript]['chrm']
    strand = data[transcript]['strand']
    protein_id = data[transcript]['protein_id']
    for cdds in data[transcript]['cds']:
        cds += cdds['seq']
    subs_dict = get_alt(cds,data[transcript])
    all_substitutions_dict[transcript]['chrm'] = chrm
    all_substitutions_dict[transcript]['strand'] = strand
    all_substitutions_dict[transcript]['protein_id'] = protein_id
    all_substitutions_dict[transcript]['cds'] = cds
    all_substitutions_dict[transcript]['all_subs'] = subs_dict


with open("all_substitutions.csv", 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['transcript', 'chrm', 'strand','protein_id','ref_nuc','coordinate','alt_nuc','aa_pos','ref_aa','alt_aa'])
    for transcript in all_substitutions_dict.keys():
        for all_sub, chan in all_substitutions_dict[transcript]['all_subs'].items():
            writer.writerow([transcript,
                            all_substitutions_dict[transcript]['chrm'],
                            all_substitutions_dict[transcript]['strand'],
                            all_substitutions_dict[transcript]['protein_id'],
                            all_sub[2],
                            chan['coordinate'],
                            all_sub[3],
                            chan['aa_pos'],
                            chan['ref_aa'],
                            chan['alt_aa']])
           





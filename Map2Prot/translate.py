from asyncio import proactor_events
import os
import sys

def get_fasta_dict(fasta_file):
    fasta_dict = {}
    
    with open(fasta_file,'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.split(">")[1].strip()
                fasta_dict[header] = ''
            else:
                fasta_dict[header] += line.strip()
            
    f.close()
    return fasta_dict

def write_fasta(fasta_dict,outfile):
    with open(outfile.split("_pan_transcipts.fasta")[0]+"_pan_proteins.fasta",'w') as f:
        for k in fasta_dict.keys():
            f.write(">" + k + "\n" + fasta_dict[k] + "\n")

    f.close()


def get_gap_positions(sequence):
    gap_positions = []
    non_gap_positions = []
    for i in range(0,len(sequence)):
        if sequence[i] == "-":
            gap_positions.append(i)
        if sequence[i] != "-":
            non_gap_positions.append(i)

    return gap_positions, non_gap_positions


def translate(transcript_file):
    transcript_dict = get_fasta_dict(transcript_file)
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        '---':'-'
    }

    protein_dict = {}
    for k,v in transcript_dict.items():
        print(k)
        protein =""
        translate_dict = {}
        if len(v)%3 == 0:
            gap_positions, non_gap_positions = get_gap_positions(v)
            print(non_gap_positions)
            print("****")
            for non_gap_pos in range(0, len(non_gap_positions), 3):
                a, b, c = non_gap_positions[non_gap_pos:non_gap_pos+3]
                codon = v[a] + v[b] + v[c] 
                codon = codon.upper()
                translate_dict[non_gap_positions[non_gap_pos]] = table[codon]
            for gap_pos in range(0, len(gap_positions), 3):
                f,s,t = gap_positions[gap_pos:gap_pos+3]
                codon = v[f] + v[s] + v[t]
                codon = codon.upper()
                translate_dict[gap_positions[gap_pos]] = table[codon]
        sorted_positions =sorted(list(translate_dict.keys()))
        for p in sorted_positions:
            protein += translate_dict[p]
                    
        protein_dict[k] = protein
    
    write_fasta(protein_dict,transcript_file)


if __name__ == "__main__":
    pan_isoforms_file = sys.argv[1]
    translate(pan_isoforms_file)



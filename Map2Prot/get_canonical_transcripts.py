from get_cds import *
from main_translate import *
from fasta_dict import *
import json

cds_dict = get_all_isoform_cds("/cta/users/abircan/pan_isoform/resources/GTFs_release_106/Homo_sapiens/Homo_sapiens.gtf", "/cta/users/abircan/pan_isoform/resources/Genomes_release_106/Homo_sapiens/Homo_sapiens.fa")


protein_to_transcript_dict = {}
with open('HUMAN_9606_idmapping.dat', 'r') as fp:
    lines = fp.readlines()
    for line in lines:
        if "ENST" in line:
            if "-" in line.split("\t")[0]:
                if  not line.split("\t")[0][:-2] in protein_to_transcript_dict.keys():
                    protein_to_transcript_dict[line.split("\t")[0][:-2]] = [line.split("\t")[-1].split(".")[0].strip()]
                else:
                    protein_to_transcript_dict[line.split("\t")[0][:-2]].append(line.split("\t")[-1].split(".")[0].strip())
            else:
                if  not line.split("\t")[0] in protein_to_transcript_dict.keys():
                    protein_to_transcript_dict[line.split("\t")[0]] = [line.split("\t")[-1].split(".")[0].strip()]
                else:
                    protein_to_transcript_dict[line.split("\t")[0]].append(line.split("\t")[-1].split(".")[0].strip())


all_proteins_dict = {}

for gene in cds_dict.keys():
    for transcript in cds_dict[gene]['transcripts'].keys():
        print(transcript)
        q_transcript = transcript
        cds = get_single_transcript_cds(cds_dict,q_transcript)
        protein_seq = translate(cds)
        all_proteins_dict[q_transcript] = {'strand':cds_dict[gene]['strand'],'chrm':cds_dict[gene]['chrm'],'cds':cds_dict[gene]['transcripts'][q_transcript],'protein_seq':protein_seq
        }
        
#print(all_proteins_dict)          
canonical_fasta_dict = get_fasta_dict("canonical_Homo_Sapiens_proteins.fasta")

canonical_transcript_dict = {}
id_not_found_dict = {}
seq_not_match_dict = {}

with open("exceptional_not_found_in_id_mapping","w") as id_f:
    for k,v in canonical_fasta_dict.items():
        for transcript in all_proteins_dict.keys():
            if all_proteins_dict[transcript]['protein_seq'] == v:
                print(v)
                if k.split("|")[1].split("|")[0] in protein_to_transcript_dict.keys():
                    if transcript in protein_to_transcript_dict[k.split("|")[1].split("|")[0]]:
                        if transcript not in canonical_transcript_dict.keys():
                            canonical_transcript_dict[transcript] = {'strand':all_proteins_dict[transcript]['strand'],'chrm':all_proteins_dict[transcript]['chrm'],'cds':all_proteins_dict[transcript]['cds'], 'protein_id':k,'protein_seq':v}
                
                else:
                    id_not_found_dict[transcript] = {'strand':all_proteins_dict[transcript]['strand'],'chrm':all_proteins_dict[transcript]['chrm'],'cds':all_proteins_dict[transcript]['cds'], 'protein_id':k,'protein_seq':v}
                    id_f.write(k + "\n")
                
for k,v in canonical_fasta_dict.items():
    for transcript in all_proteins_dict.keys():
        if all_proteins_dict[transcript]['protein_seq'] != v:
            if transcript not in id_not_found_dict.keys() or transcript not in canonical_transcript_dict.keys() or transcript not in seq_not_match_dict.keys():
                seq_not_match_dict[transcript] =  {'strand':all_proteins_dict[transcript]['strand'],'chrm':all_proteins_dict[transcript]['chrm'],'cds':all_proteins_dict[transcript]['cds'], 'protein_id':k,'protein_seq':v}


with open('all_canonical_transcripts.json', 'w') as fp:
    json.dump(canonical_transcript_dict, fp)


with open('all_id_not_found.json', 'w') as fp1:
    json.dump(id_not_found_dict, fp1)


with open('exceptional_sequence_not_matched.json', 'w') as fp2:
    json.dump(seq_not_match_dict, fp2)

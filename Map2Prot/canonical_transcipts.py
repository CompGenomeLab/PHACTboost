from get_cds import *
from main_translate import *
from fasta_dict import *
import json


def get_transcript_from_protein(id_mapping_file):
    protein_to_transcript_dict = {}
    with open(id_mapping_file, 'r') as fp:
        lines = fp.readlines()
        for line in lines:
            if "Ensembl_TRS" in line:
                if "-" in line.split("\t")[0]:
                    if line.split("\t")[0][-2:] == "-1":
                        if  not line.split("\t")[0][:-2] in protein_to_transcript_dict.keys():
                            protein_to_transcript_dict[line.split("\t")[0][:-2]] = [line.split("\t")[-1].split(".")[0].strip()]
                        else:
                            protein_to_transcript_dict[line.split("\t")[0][:-2]].append(line.split("\t")[-1].split(".")[0].strip())
                else:
                    if  not line.split("\t")[0][:-2] in protein_to_transcript_dict.keys():
                            protein_to_transcript_dict[line.split("\t")[0][:-2]] = [line.split("\t")[-1].split(".")[0].strip()]
                    else:
                            protein_to_transcript_dict[line.split("\t")[0][:-2]].append(line.split("\t")[-1].split(".")[0].strip())
    return protein_to_transcript_dict


def get_all_proteins(cds_dict):

    all_proteins_dict = {}

    for gene in cds_dict.keys():
        for transcript in cds_dict[gene]['transcripts'].keys():
            print(transcript)
            q_transcript = transcript
            cds = get_single_transcript_cds(cds_dict,q_transcript)
            protein_seq = translate(cds)
            all_proteins_dict[q_transcript] = {'strand':cds_dict[gene]['strand'],'chrm':cds_dict[gene]['chrm'],'cds':cds_dict[gene]['transcripts'][q_transcript],'protein_seq':protein_seq
            }
    return all_proteins_dict
        
          

def get_canonicals(all_proteins_dict,protein_to_transcript_dict,canonical_fasta_dict):
    canonical_transcript_dict = {}

    with open("not_found_in_id_mapping","w") as id_f:
        for k,v in canonical_fasta_dict.items():
            for transcript in all_proteins_dict.keys():
                if all_proteins_dict[transcript]['protein_seq'] == v:
                    if k.split("|")[1].split("|")[0] in protein_to_transcript_dict.keys():
                        if transcript in protein_to_transcript_dict[k.split("|")[1].split("|")[0]]:
                            canonical_transcript_dict[transcript] = {'strand':all_proteins_dict[transcript]['strand'],'chrm':all_proteins_dict[transcript]['chrm'],'cds':all_proteins_dict[transcript]['cds'], 'protein_id':k,'protein_seq':v}
                    else:
                        id_f.write(k + "\n")


    with open('canonical_transcripts.json', 'w') as fp:
        json.dump(canonical_transcript_dict, fp)



if __name__ == "__main__":
    gtf_file = sys.argv[1]
    genome_file = sys.argv[2]
    id_mapping_file = sys.argv[3]
    canonical_proteins_fasta = sys.argv[4]
    cds_dict = get_all_isoform_cds(gtf_file,genome_file)
    protein_to_transcript_dict = get_transcript_from_protein(id_mapping_file) 
    all_proteins_dict = get_all_proteins(cds_dict)
    canonical_fasta_dict = get_fasta_dict(canonical_proteins_fasta)
    get_canonicals(all_proteins_dict,protein_to_transcript_dict,canonical_fasta_dict)




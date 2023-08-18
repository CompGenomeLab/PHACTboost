import os
import sys
import json
from main_translate import *
import re
from math import *

complement_dict = {"A":"T","T":"A","G":"C","C":"G"}

def get_nuc(ChrPositions_file):
    chr_pos_dict = {}
    with open(ChrPositions_file,'r') as f:
        lines = f.readlines()
        for line in lines:
            if  line and "chrs" not in line:
                #print(line) 
                l = None
                if ChrPositions_file.split(".")[-1] == "csv":
                    l =  line.split(",")[:]
                elif ChrPositions_file.split(".")[-1] == "tsv":
                    l =  line.split("\t")[:]
                elif  ChrPositions_file.split(".")[-1] == "txt":
                    l =  line.split(" ")[:]
                chrm = l[0].strip()
                position = l[1].strip()
                ref_nuc = l[2].strip()
                alt_nuc = l[3].strip()
                if chrm in chr_pos_dict.keys():
                    chr_pos_dict[chrm].append({'pos':position,'ref_nuc':ref_nuc,'alt_nuc':alt_nuc})
                else:
                    chr_pos_dict[chrm] = [{'pos':position,'ref_nuc':ref_nuc,'alt_nuc':alt_nuc}]
        

    return chr_pos_dict

def get_json_dict(canonical_transcripts_json_file):
    with open(canonical_transcripts_json_file) as json_file:
        canonical_transcripts_dict = json.load(json_file)
    return canonical_transcripts_dict



def get_alt_nuc(cds,position,alt_nuc,strand,ref_nuc):
    alt_nuc_pos = 0
    all_coordinates = [c['coordinates'] for c in cds]
    cds_list = [c['seq'] for c in cds]
    original_cds = "".join([c['seq'] for c in cds]).upper()
    _ = list(original_cds)
    if strand == "-":
        for i in range(len(all_coordinates)):
            if int(all_coordinates[i][0]) <= int(position) and int(all_coordinates[i][1]) >= int(position):
                alt_nuc_pos +=  abs(int(all_coordinates[i][1])- int(position)) + 1
                break
            else:
                alt_nuc_pos += abs(int(all_coordinates[i][1]) - int(all_coordinates[i][0])) + 1 
        _[alt_nuc_pos-1] = complement_dict[alt_nuc]
        alt_cds = "".join(_)
        alt_cds = alt_cds.upper()

    elif strand == "+":
        for i in range(len(all_coordinates)):
            if int(all_coordinates[i][0]) <= int(position) and int(all_coordinates[i][1]) >= int(position):
                alt_nuc_pos +=  abs(int(all_coordinates[i][0])- int(position)) + 1
                break
            else:
                alt_nuc_pos += abs(int(all_coordinates[i][1]) - int(all_coordinates[i][0])) + 1 

        _[alt_nuc_pos-1] = alt_nuc
        alt_cds = "".join(_)
        alt_cds = alt_cds.upper()
        
    alt_aa_pos = ceil(alt_nuc_pos/3)
    if strand == "-":
        if complement_dict[original_cds[alt_nuc_pos-1]] == ref_nuc:
            return alt_nuc_pos, alt_cds, alt_aa_pos
    elif strand == "+":
        if original_cds[alt_nuc_pos-1] == ref_nuc:
            return alt_nuc_pos, alt_cds, alt_aa_pos


def get_alt(canonical_transcripts_dict,chrm,position,ref_nuc,alt_nuc):
    original_cds = ""
    original_protein_id = ""
    original_protein_seq = ""
    alt_cds = ""
    for transcript in canonical_transcripts_dict.keys():
        for values in canonical_transcripts_dict[transcript]:
            if values['chrm'] == chrm:
                for coordinates in values['cds']:
                    coord =  [int(i) for i in coordinates['coordinates']]
                    if int(position) >= coord[0] and  coord[1] >= int(position):
                        strand = values['strand']
                        alt_nuc_pos, alt_cds, alt_aa_pos = get_alt_nuc(values['cds'],position,alt_nuc,strand,ref_nuc)
                        all_cds = [cds['seq'] for cds in values['cds']]
                        original_cds = "".join(all_cds)
                        original_protein_id = values['protein_id']
                        original_protein_seq = values['protein_seq']

                           
    original_cds = original_cds.upper()
    if original_protein_id:
        return original_cds,alt_cds, original_protein_id, original_protein_seq, alt_nuc_pos, alt_aa_pos
   

                     
def get_protein_pos(chr_pos_dict,canonical_transcripts_dict,mapped_file,not_found_file):
    not_match = {}
    with open(mapped_file,'w') as f:
        f.write("chrm" + "\t" + "position" + "\t" + "cds_position" + "\t" + "ref_nuc"+ "\t" +"alt_nuc" + "\t"  +"alt_protein_position"+ "\t" + "original_aa" + "\t"+ "alt_aa"+ "\t" +"original_protein_id"+ "\t" +"mutation_type"+ "\n")
        for chrm, v in chr_pos_dict.items():
            for alternation in v:
                pos = alternation['pos']
                ref_nuc = alternation['ref_nuc']
                alt_nuc = alternation['alt_nuc']
                try:
                    original_cds,alt_cds, original_protein_id, original_protein_seq, alt_nuc_pos, alt_aa_pos =  get_alt(canonical_transcripts_dict,chrm,int(pos),ref_nuc,alt_nuc)
                    alt_protein_seq = translate_substitution(alt_cds)
                    alt_aa = alt_protein_seq[alt_aa_pos-1]
                    org_aa = original_protein_seq[alt_aa_pos-1]
                    if "*" in alt_aa:
                        f.write(chrm + "\t" + pos + "\t" + str(alt_nuc_pos) +  "\t" + ref_nuc+ "\t" + alt_nuc +  "\t"   + str(alt_aa_pos) + "\t" +org_aa + "\t" + alt_aa+  "\t" + original_protein_id +"\t" + "early_stop_codon"+"\n")                 

                    elif len(original_cds) < len(alt_cds):
                        f.write(chrm + "\t" + pos + "\t" + str(alt_nuc_pos) +  "\t" + ref_nuc+ "\t" + alt_nuc +  "\t"   + str(alt_aa_pos) + "\t" +org_aa + "\t" + alt_aa+  "\t" + original_protein_id +"\t" + "insertion"+"\n")                 


                    elif len(original_cds) > len(alt_cds):
                        f.write(chrm + "\t" + pos + "\t" + str(alt_nuc_pos) +  "\t" + ref_nuc+ "\t" + alt_nuc +  "\t"   + str(alt_aa_pos) + "\t" +org_aa + "\t" + alt_aa+  "\t" + original_protein_id +"\t" + "deletion"+"\n")                 

                    elif  alt_aa == org_aa:

                        f.write(chrm + "\t" + pos + "\t" + str(alt_nuc_pos) +  "\t" + ref_nuc+ "\t" + alt_nuc +  "\t"   + str(alt_aa_pos) + "\t" +org_aa + "\t" + alt_aa+  "\t" + original_protein_id +"\t" + "synonymous"+"\n")                 

                    else:
                        f.write(chrm + "\t" + pos + "\t" + str(alt_nuc_pos) +  "\t" + ref_nuc+ "\t" + alt_nuc +  "\t"   + str(alt_aa_pos) + "\t" +org_aa + "\t" + alt_aa+  "\t" + original_protein_id +"\t" + "non_synonymous"+"\n")                 
                except: 
                    not_match[chrm] = v
                    
    f.close()
    with open(not_found_file,"w") as file:
        file.write("chrm" + "\t" + "position" + "\t" + "ref_nuc"+ "\t" +"alt_nuc"+ "\n")
        for chrm, v in chr_pos_dict.items():
            for alternation in v:
                pos = alternation['pos']
                ref_nuc = alternation['ref_nuc']
                alt_nuc = alternation['alt_nuc']
                file.write(chrm + "\t" + ref_nuc + "\t" + alt_nuc+"\n")
    file.close()
                            



if __name__ == "__main__":
    variants_file = sys.argv[1]
    canonical_transcripts_json = sys.argv[2]
    mapped_file = sys.argv[3]
    not_found_file = sys.argv[4]
    chr_pos_dict= get_nuc(variants_file)
    canonical_transcripts_dict = get_json_dict(canonical_transcripts_json)
    get_protein_pos(chr_pos_dict,canonical_transcripts_dict,mapped_file,not_found_file)





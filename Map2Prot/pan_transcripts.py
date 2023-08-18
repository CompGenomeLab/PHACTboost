import subprocess
import os
from pprint import pprint
import sys
import subprocess
import os
from pprint import pprint
import sys


def get_all_isoform_cds(gtf_file,genome_file):
    cds_dict = {}
    transcript_types = ["protein_coding","IG_V_gene", "TR_V_gene", "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", "TR_C_gene", "TR_J_gene", "TR_D_gene","ambiguous_orf","translated_processed_pseudogene","translated_unprocessed_pseudogene","TEC"]

    with open(gtf_file,'r') as f:
        for line in f:
            if not line.startswith("#"):
                fields = line.rstrip().split('\t')
                attributes = fields[8].replace('"',"")
                gene_id = attributes.split("gene_id ")[1].split(";")[0]
                type = fields[2]
                chrm = fields[0]
                strand = fields[6]
                if type == "CDS":
                    exon =  attributes.split("exon_number ")[1].split(";")[0]
                    start = fields[3]
                    stop = fields[4]
                    transcript_id = attributes.split("transcript_id ")[1].split(";")[0]
                    transcript_bio_type = attributes.split("transcript_biotype ")[1].split(";")[0]
                    if transcript_bio_type in transcript_types:
                        if strand=="-":                            
                            output = os.popen("samtools faidx " + genome_file + " " +chrm + ":" + start + "-" + stop + " -i").read().split("\n")
                            seq = "".join(output[1:])
                        if strand=="+":
                            output = os.popen("samtools faidx " + genome_file + " " +chrm + ":" + start + "-" + stop).read().split("\n")
                            seq = "".join(output[1:])
                        if not cds_dict:
                            cds_dict[gene_id] = {'chrm':chrm,'strand':strand,'transcripts':{transcript_id:[{'exon':exon,'coordinates':(start,stop),'seq':seq}]}}
                        else:
                            if gene_id in cds_dict.keys():
                                if transcript_id not in cds_dict[gene_id]['transcripts'].keys():
                                    cds_dict[gene_id]['transcripts'].update({transcript_id:[{'exon':exon,'coordinates':(start,stop),'seq':seq}]})
                                else:
                                    cds_dict[gene_id]['transcripts'][transcript_id].append({'exon':exon,'coordinates':(start,stop),'seq':seq})
                            else:
                                cds_dict[gene_id] = {'chrm':chrm,'strand':strand,'transcripts':{transcript_id:[{'exon':exon,'coordinates':(start,stop),'seq':seq}]}}

    return cds_dict


#cds_dict=get_all_isoform_cds("/cta/users/abircan/pan_isoform/resources/GTFs_release_106/Homo_sapiens/Homo_sapiens.gtf","/cta/users/abircan/pan_isoform/resources/Genomes_release_106/Homo_sapiens/Homo_sapiens.fa")
def get_single_transcript_cds(cds_dict,transcript_id):
    cds = ""
    for gene in cds_dict.keys():
        for transcript in cds_dict[gene]['transcripts'].keys():
            if transcript == transcript_id:
                for exons in cds_dict[gene]['transcripts'][transcript_id]:
                    cds +=  exons['seq']
    return cds.upper()


def get_coordinates(cds_dict):
    coordinates_dict = {}
    for gene in cds_dict.keys():
        coordinates = []
        strand = cds_dict[gene]['strand']
        for transcript,exons in cds_dict[gene]['transcripts'].items():
            coordinates += ([exon['coordinates'] for exon in exons])
        coordinates = list(dict.fromkeys(coordinates))
        if strand=="+":
            coordinates = sorted(coordinates, key=lambda x: x[0])
        elif strand=="-":
            coordinates = sorted(coordinates, key=lambda x: x[0],reverse=True)    
        coordinates_dict[gene] = coordinates
    return coordinates_dict

def get_pan_transcript_seqs(cds_dict,coordinates_dict):
    pan_transcript_dict = {}
    for gene in cds_dict.keys():
        pan_transcript_dict[gene] = {}
        coordinates = coordinates_dict[gene]
        for transcript,exons in cds_dict[gene]['transcripts'].items():
            pan_transcript_dict[gene][transcript] = ""
            for coordinate in coordinates:
                if any(coordinate==exon['coordinates'] for exon in exons):
                    for exon in exons:
                        if coordinate == exon['coordinates']:
                            pan_transcript_dict[gene][transcript] += exon['seq']
                else:
                    pan_transcript_dict[gene][transcript] += "-"*(abs((int(coordinate[1])-int(coordinate[0])))+1)
    return pan_transcript_dict

def write_fasta(pan_transcript_dict,outfile):
    with open(outfile,'w') as f:
        for gene in pan_transcript_dict.keys():
            for transcript in pan_transcript_dict[gene].keys():
                f.write(">" + "gene_" + gene + "_" + "transcript_"+ transcript + "\n" + pan_transcript_dict[gene][transcript]+"\n")
    f.close()


# # if __name__ == "__main__":
# #     gtf_file = sys.argv[1]
# #     genome_file = sys.argv[2]
# #     outfile = sys.argv[3]
# #     cds_dict = get_all_isoform_cds(gtf_file,genome_file)
# #     coordinates_dict = get_coordinates(cds_dict)
# #     pan_transcript_dict = get_pan_transcript_seqs(cds_dict,coordinates_dict)
# #     write_fasta(pan_transcript_dict,outfile)



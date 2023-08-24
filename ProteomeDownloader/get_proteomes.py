import urllib.request
import gzip
import shutil
import csv
import time
start_time = time.time()
#Proteome ID file accesssion date and url
#time of download 19-07-2021 14:29


#provide raw URL, XXX will be replaced with proteome ID, YYY will be replaced by organism ID
raw_url = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/XXX/XXX_YYY.fasta.gz'
#read the tab delimited file containing proteome IDs and organism IDs
with open("/cta/users/abircan/check/proteomes_list.tab", newline="") as file:
    file_reader = csv.DictReader(file, delimiter = "\t")
    for line in file_reader:
    	#get required identifiers
        proteome_id = line["Proteome ID"]
        organism = line["Organism"].replace(" ", "_") #remove gaps
        organism = organism.replace("/", "-") #/ sign is considered as directory
        organism = organism[:24] #too long names cause errors
        organism_id = line["Organism ID"]
        #provide url to access reference proteome
        currenturl = raw_url.replace("XXX", proteome_id)
        currenturl = currenturl.replace("YYY", organism_id)
        #provide fasta file name (includes important identifiers)
        fastaname = proteome_id + "_" + organism_id + "_" + organism + ".fasta"
        #provide path to download compresssed files
        downpath = "/cta/users/abircan/check/database/compressed_files/"
        #request the proteome, download it to the provided path
        try:
            urllib.request.urlretrieve(currenturl,downpath + fastaname + ".gz")
            #provide path to unzip files
            unzippath = "/cta/users/abircan/check/database/canonical_reference_proteomes/"
            #open the file, unzip it, write it down, close the file
            with gzip.open(downpath + fastaname + ".gz", 'rb') as f_in:
                with open(unzippath + fastaname, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        except:
            print(fastaname)
            continue
print("My program took", time.time() - start_time, "seconds to run")


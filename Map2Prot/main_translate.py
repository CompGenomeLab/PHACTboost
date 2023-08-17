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
    #print(s[:length_divisible_by_three])
    return s[:length_divisible_by_three]

#longest_substring_divisible_by_3("abcde")

def translate(ORF):
    ORF = longest_substring_divisible_by_3(ORF)
    aa = ''
    for i in range(0,len(ORF),3):
        dictkey=ORF[i:i+3].replace("T","U")
        aa+=universal[dictkey]
    if aa:
        aa= aa[:-1]
        _ = list(aa)
        _[0] = "M"
        #print("".join(_))
        return "".join(_)


def translate_substitution(ORF):
    ORF = longest_substring_divisible_by_3(ORF)
    aa = ''
    for i in range(0,len(ORF),3):
        dictkey=ORF[i:i+3].replace("T","U")
        aa+=universal[dictkey]
    if aa:
        aa= aa[:-1]
        return aa



#longest_substring_divisible_by_3("ATGCTGCTGCTTCTGCTGCTTCTGGGGCCAGGCTCCGGGCTTGGTGCTGTCGTCTCTCAACATCCGAGCAGGGTTATCTGTAAGAGTGGAACCTCTGTGAAGATCGAGTGCCGTTCCCTGGACTTTCAGGCCACAACTATGTTTTGGTATCGTCAGTTCCCGAAACAGAGTCTCATGCTGATGGCAACTTCCAATGAGGGCTCCAAGGCCACATACGAGCAAGGCGTCGAGAAGGACAAGTTTCTCATCAACCATGCAAGCCTGACCTTGTCCACTCTGACAGTGACCAGTGCCCATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTAGAGA")

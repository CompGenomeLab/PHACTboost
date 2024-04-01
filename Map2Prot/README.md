# Genomic and Proteomic Analysis Toolkit

## Overview

These Python scripts extract coding sequences from genomic data, find their genomic coordinates, exon, strand, chromosome information and translate these sequences to proteins, match them with canonical protein sequences. All possible single nucleotide substitutions are calculated and corresponding amino acids and protein positions are found. 


- **fasta_dict.py**: Generates a dictionary from a FASTA file, associating sequence headers with their corresponding sequences.

- **main_translate.py**: Translates nucleotide sequences to amino acid sequences.

- **get_cds.py**: Extracts CDS from GTF and genome FASTA files, producing a dictionary mapping gene IDs to their CDS information.

- **get_canonical_transcripts.py**: Matches translated protein sequences with canonical protein sequences, outputting a JSON file with transcripts whose translated proteins align with canonical sequences. Additional JSON files for IDs not found and non-matching sequences are also generated.

- **prot2transcript.py**: Calculates all potential single nucleotide substitutions for canonical transcripts, outputting a CSV file including genomic coordinates, amino acid positions, reference nucleotide, alternatig nucleotide, reference amino acid, altered amino acid etc

## Workflow

1. **CDS Extraction**: Use `get_cds.py` with GTF and genome FASTA files to extract CDS information.

2. **Canonical Transcripts Identification**: Execute `get_canonical_transcripts.py`, which requires the CDS output, a canonical human proteins FASTA file, and an ID mapping file to match protein sequences with transcripts.

3. **Amino Acid Substitutions Analysis**: Run `prot2transcript.py` to identify all possible single nucleotide substitutions for canonical transcripts, generating a detailed CSV file.

## Requirements

- Python 3.x.
- External libraries: `pandas` (for `prot2transcript.py`).
- Samtools (for sequence extraction in `get_cds.py`).
- Ensure all prerequisite files (e.g., GTF, genome FASTA, canonical proteins FASTA, ID mapping files) are accessible.


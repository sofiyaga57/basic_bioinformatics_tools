# bio_files_processor.py

**basic_bioinformatic_tools** is a tool which allows fresh 
bioinformatics to work with fasta, fastq and gbk files. 

This script has three main functions:
1) `convert_multiline_fasta_to_oneline`
Converts multiline fasta file to oneline fasta file. 
Returns file, names of sequences are seq1, seq2, etc.
2) `select_genes_from_gbk_to_fasta`
Selects genes surrounding genes of interest in the gbk file and return their sequences in 
a new fasta file.
3) `change_fasta_start_pos`


# basic_bioinformatics_tools.py

**basic_bioinformatic_tools** is a tool which allows fresh bioinformatics to try their hand in working with DNA, RNA, protein sequences, as well as with fastq files-dictionaries.

This tools consists of three modules discussed in great detail below:
1) `dna_rna_tools.py`
2) `protein_tools.py`
3) `fastq_tools.py`

Script `basic_bioinformatics_tools.py` imports three main functions from these three modules.

## dna_rna_tools.py

### Usage

This utility is realized inside the file `dna_rna_tools.py`. This program contains the main function `run_dna_rna_tools`, as well as other additional functions described below. 
The function takes as input an arbitrary number of arguments with DNA or RNA sequences in str format, the last argument being the name of the str-formatted procedure to be implemented.
If one sequence is input, a string with the result is returned, if several sequences are input, an array of strings with the result is returned. 
If a string that does not match DNA or RNA is input (the sequence contains letters other than ATGCU, or the sequence contains U and T at the same time), ValueError: Not a nucleic acid is returned. 
The utility does not change the case of the characters in the sequence.

### Procedures

- `transcribe` - print the transcribed DNA. Outputs an RNA sequence. Gives ValueError: Not a DNA if 
the sequence is not a DNA.
- `reverse_transcribe` - print the reverse transcribed RNA. Same as transcribe, but makes DNA from RNA. 
Raises ValueError if not an RNA is given.
- `reverse` - print the reversed sequence. Raises ValueError: Not a DNA or RNA if not a DNA or RNA is input.
- `complement` - print complementary sequence. Raises ValueError: Not a DNA or RNA if not a DNA or RNA is input.
- `reverse_complement` - print reverse complementary sequence, same as reverse+complement
- `compute_melting_temperature` - print DNA melting temperature value, works only with DNA sequences, otherwise 
raises ValueError: Not a DNA.

### Examples

```python
run_dna_rna_tools('ATG', 'transcribe') 
# 'AUG'
run_dna_rna_tools('ATG', 'TTg', 'transcribe') 
# ['AUG', 'UUg']
run_dna_rna_tools('AUG', 'reverse_transcribe') 
# ['ATG']
run_dna_rna_tools('ATg', 'reverse') 
# 'gTA'
run_dna_rna_tools('ATG', 'TTg', 'complement') 
# ['TAC', 'AAc']
run_dna_rna_tools('ATg', 'reverse_complement') 
# 'cAT'
run_dna_rna_tools('ATGATG', 'compute_melting_temperature') 
# 20.0
```

## protein_tools.py

**protein_tools.py** - is a tool which allows the performing of various procedures for a user entered protein sequences. 

### Usage

The tool works by calling the function `run_protein_tools`, which takes arbitrary number of arguments with protein 
sequencies (*str*) and the name of the procedure to be performed (always the last argument, *str*, see the usage 
examples below). The output is the result of the procedure as *string, tuple* or *dictionary* if one sequence is 
submitted or *list* if several.

**NOTE:**  For the procedure `check_mutations` a fixed number of string arguments are used: one RNA sequence, one 
protein sequence and the name of procedure itself.

### Procedures

- `compute_molecular_weight` — computes molecular weight of protein sequence in g/mol
- `compute_length` — computes the number of amino acids in protein sequence
- `compute_hydrophobicity` — computes the percentage of hydrophobic aminoacids in protein sequence
- `check_mutations` — checks missense mutations in the protein sequence after translation
- `protein_to_dna`- returns possible variants of DNAs for a given protein sequence
- `get_frequencies_of_amino_acids` - calculates the number of each aminoacid in protein sequence

### Examples
```python
run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'compute_length')
#[10, 18]

run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'compute_molecular_weight')
#[1055.496, 1886.872]

run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'compute_hydrophobicity')
#[50.0, 27.778]

run_protein_tools('AUGGAUCAUcAAUAA', 'MDKL*', 'check_mutations')
#'Mutations: K3, L4.'

run_protein_tools('MAEGLP', 'LYGSQT','protein_to_dna')
#['ATG GCT/GCC/GCA/GCG GAA/GAG GGT/GGC/GGA/GGG TTA/TTG/CTT/CTC/CTA/CTG CCT/CCC/CCA/CCG',
#'TTA/TTG/CTT/CTC/CTA/CTG TAT/TAC GGT/GGC/GGA/GGG TCT/TCC/TCA/TCG/AGT/AGC CAA/CAG ACT/ACC/ACA/ACG']

run_protein_tools('MAEGLP', 'LYGSQT','get_frequencies_of_amino_acids')
#[{'M': 1, 'A': 1, 'E': 1, 'G': 1, 'L': 1, 'P': 1},
#{'L': 1, 'Y': 1, 'G': 1, 'S': 1, 'Q': 1, 'T': 1}]
```
   
### Additional information
- The program works **only** with protein and RNA sequences. If any of the entered sequences contain inappropriate 
characters or cannot exist, the program will display an error. Sequences can contain characters of any case.

```python
run_protein_tools('PROTEIN', 'compute_molecular_weight') #ValueError: Invalid protein sequence
run_protein_tools('AUGGAU_AUcAAUAA', 'MDKL*', 'check_mutations') #ValueError: Invalid RNA sequence
```
- For the procedure `check_mutations` there are extra requirements for RNA and protein sequences: mRNA sequences must 
contain **start-codon** and **one of the stop-codons**, protein sequnces must start with **"M"** and ends with **"*"** (stop-codon). 
```python
run_protein_tools("AUGGUAGGGAAAUUUUGA", "MGGKF", 'check_mutations') #ValueError: Stop (*) is absent
run_protein_tools("AUGGUAGGGAAAUUUUGA", "GGKF*", 'check_mutations') #ValueError: Start (M) is absent
```

## fastq_tools.py

### Usage
This tool enables one to work with fastq files and filtrate them by given conditions.
The main function run_fastq_tools takes fastq file as first argument, and output_filename as 
a second argument, the following three arguments has default meaning:
1) `gc_bounds`, sets the upper limit for filtration if only one number is given,
otherwise sets range for filtration. Default meaning is None.
2) `length_bounds`, same as gc_bounds, but default value is (0, 2**32).
3) `quality_threshold`, lower limit for filtration. Default value is 0. 

Raises ValueError("Too strict conditions") if the return dictionary in any of the functions is empty.

### Procedures

- `filter_gc_content`, filters fastq dictionary by GC-content
- `filter_length`, filters fastq dictionary by length
- `filter_quality`, filters fastq dictionary by quality

*Author contributions:* <br>
Sofiya Vinogradova: module `dna_rna_tools.py`, module `fastq_tools.py`, module `protein_tools.py` – functions `compute_length`, `get_frequencies_of_amino_acids`, `protein_to_dna` <br> 
Nikita Zherko: module `protein_tools.py` – functions `compute_hydrophobicity`, `translate_rna`, `check_mutations`
Artyom Toropov: module `protein_tools.py` – functions `is_protein`, `is_rna`, `compute_molecular_weight`, `run_protein_tools` <br> 

### Contacts
Please use contacts below to reach out with any comments, concerns, or discussions regarding **protein_tools.py.** <br>
- Sofiya Vinogradova ([@sofiyaga57](https://github.com/sofiyaga57/)) <br>

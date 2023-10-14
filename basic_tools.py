"""
This script imports three modules and their main functions.
"""

from basic_bioinformatics_tools.dna_rna_tools import run_dna_rna_tools
# from basic_bioinformatics_tools.protein_tools import run_protein_tools as run_protein_tools
# from basic_bioinformatics_tools.fastq_tools import run_fastq_tools as run_fastq_tools

print(run_dna_rna_tools('atg', 'reverse'))
# run_protein_tools()
# run_fastq_tools()

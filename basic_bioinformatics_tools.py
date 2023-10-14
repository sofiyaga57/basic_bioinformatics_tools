"""
This script imports three modules and their main functions.
"""

import basic_bioinformatics_tools.dna_rna_tools as dna_rna_tools
import basic_bioinformatics_tools.protein_tools as protein_tools
import basic_bioinformatics_tools.fastq_tools as fastq_tools

dna_rna_tools.run_dna_rna_tools()
protein_tools.run_protein_tools()
fastq_tools.run_fastq_tools()

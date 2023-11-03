DNA_NUCLEOTIDES = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}


RNA_NUCLEOTIDES = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}


def is_dna(seq: str) -> bool:
    """
    Checks whether sequence is DNA.
    :param seq: str, nucleotide sequence
    :return: bool
    """
    unique_chars = set(seq)
    nucleotides = set('ATGCatgc')
    return unique_chars.issubset(nucleotides)


def is_rna(seq: str) -> bool:
    """
    Checks whether sequence is RNA.
    :param seq: str, nucleotide sequence
    :return: bool
    """
    unique_chars = set(seq)
    nucleotides = set('AUGCaugc')
    return unique_chars.issubset(nucleotides)


def transcribe(*seqs: str) -> list:
    """
    Transcribes DNA into RNA. Changes T/t to U/u nucleotides.
    :param seqs: str, DNA sequences
    :return: str if one sequence, else returns list of sequences
    Raises ValueError('Not a DNA') if sequence is not a DNA
    """
    transcribed_seqs = []
    for seq in seqs:
        if not is_dna(seq):
            raise ValueError('Not a DNA')
        transcribed_seqs.append(seq.replace('T', 'U').replace('t', 'u'))
    return transcribed_seqs if len(transcribed_seqs) > 1 else transcribed_seqs[0]


def reverse_transcribe(*seqs: str) -> list:
    """
    Transcribes RNA into DNA. Changes U/u to T/t nucleotides.
    :param seqs: str, RNA sequences
    :return: str if one sequence is given, else returns list of sequences
    Raises ValueError('Not a DNA') if sequence is not a RNA
    """
    reverse_transcribed_seqs = []
    for seq in seqs:
        if not is_rna(seq):
            raise ValueError('Not a RNA')
        reverse_transcribed_seqs.append(seq.replace('U', 'T').replace('u', 't'))
    return reverse_transcribed_seqs if len(reverse_transcribed_seqs) > 1 else reverse_transcribed_seqs[0]


def reverse(*seqs):
    """
    Return reverse sequence to a given DNA or RNA sequence.
    :param seqs: str, DNA or RNA sequence
    :return: str if one sequence is given, else returns list of sequences
    Raises ValueError('Not a DNA or RNA') if sequence is not a DNA or RNA
    """
    reverse_seqs = []
    for seq in seqs:
        if not is_dna(seq) and not is_rna(seq):
            raise ValueError('Not a DNA or RNA')
        reverse_seqs.append(''.join(list(seq)[::-1]))
    return reverse_seqs if len(reverse_seqs) > 1 else reverse_seqs[0]


def complement(*seqs: str) -> list:
    """
    Returns complement to the input DNA or RNA sequence
    :param seqs: str, DNA or RNA sequence
    :return: str if one sequence is given, else returns list of sequences
    Raises ValueError('Not a DNA or RNA') if sequence is not a DNA or RNA
    """
    complement_seqs = []
    for seq in seqs:
        if is_dna(seq):
            complement_dna = ''
            for nucleotide in seq:
                complement_dna += DNA_NUCLEOTIDES[nucleotide]
            complement_seqs.append(complement_dna)
        elif is_rna(seq):
            complement_rna = ''
            for nucleotide in seq:
                complement_rna += RNA_NUCLEOTIDES[nucleotide]
            complement_seqs.append(complement_rna)
        else:
            raise ValueError('Not a DNA or RNA')
    return complement_seqs if len(complement_seqs) > 1 else complement_seqs[0]


def reverse_complement(*seqs: str) -> list:
    """
    Returns reverse complement to the input DNA or RNA sequence
    :param seqs: str, DNA or RNA sequence
    :return: str if one sequence is given, else returns list of sequences
    Raises ValueError('Not a DNA or RNA') if sequence is not a DNA or RNA
    """
    reverse_complement_seqs = []
    for seq in seqs:
        if is_dna(seq) or is_rna(seq):
            reverse_complement_seqs.append(''.join(complement(reverse(seq))))
        else:
            raise ValueError('Not a DNA or RNA')
    return reverse_complement_seqs if len(reverse_complement_seqs) > 1 else reverse_complement_seqs[0]


def compute_melting_temperature(*seqs: str) -> float:
    """
    Computes melting temperature of the input DNA sequence
    :param seqs: str, DNA sequence
    :return: float, melting temperature, rounded to 1 decimal place
    Raises ValueError('Not a DNA') if sequence is not a DNA.
    Raises ValueError('Wrong length') if sequence's length is less than 6 and more than 50 nucleotides.
    """
    melting_temperatures = []
    for seq in seqs:
        if not is_dna(seq):
            raise ValueError('Not a DNA')
        if 6 <= len(seq) < 14:
            melting_temperatures.append(round(float((seq.upper().count('A') + seq.upper().count('T')) * 2
                                                    + (seq.upper().count('G') + seq.upper().count('C')) * 3),
                                              ndigits=1))
        elif 14 <= len(seq) <= 50:
            melting_temperatures.append(round(64.9 + 41 * float((seq.upper().count('G') +
                                                                 seq.upper().count('C') - 16.4)) / len(seq),
                                              ndigits=1))
        else:
            raise ValueError('Wrong length')
    return melting_temperatures if len(melting_temperatures) > 1 else melting_temperatures[0]


def run_dna_rna_tools(*args: str) -> list:
    """
    Runs dna_rna_tools, available functions: transcribe, reverse_transcribe, reverse, reverse_complement,
    compute_melting_temperature.
    :param args: str, DNA or RNA sequences, last argument stand for desired function
    :return: str if one sequence is given, else returns list of sequences
    Raises ValueErrors depending on the function applied.
    """
    list_of_functions = {'reverse': reverse, 'complement': complement, 'transcribe': transcribe,
                         'reverse_transcribe': reverse_transcribe, 'reverse_complement': reverse_complement,
                         'compute_melting_temperature': compute_melting_temperature}
    *seqs, function = args
    return list_of_functions[function](*seqs)

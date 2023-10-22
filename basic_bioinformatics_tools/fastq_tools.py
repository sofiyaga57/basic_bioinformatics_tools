import os


def compute_gc_content(seq: str) -> float:
    """
    Computes GC-content of the input sequence.
    :param seq: str, nucleotide sequence
    :return: float, result of computation
    """
    return ((seq.count('G') + seq.count('C')) / len(seq)) * 100


def compute_nucleotide_quality(seq: str) -> float:
    """
    Computes average nucleotide phred33 quality.
    :param seq: str, quality sequence
    :return: int, computed average quality
    """
    quality = 0
    for nucleotide in list(seq):
        quality += ord(nucleotide) - 33
    return quality/len(seq)


def filter_gc_content(seqs: dict, gc_bounds=None) -> dict:
    """
    Filters fastq dictionary by GC-content.
    :param seqs: dict, fastq dictionary.
    :param gc_bounds: float if one value is given (upper limit of filtration),
    tuple – otherwise (bounds of filtration). If no arguments are given, default value is None, returns input dictionary.
    :return: dict, filtered fastq dictionary.
    Raises ValueError("Too strict conditions") if the return dictionary is empty.
    """
    if gc_bounds is None:
        return seqs
    seqs_filtered: dict = {}
    for name, seq in seqs.items():
        if isinstance(gc_bounds, tuple):
            if gc_bounds[0] <= compute_gc_content(seq[0]) <= gc_bounds[1]:
                seqs_filtered[name] = seq
        else:
            if compute_gc_content(seq[0]) <= gc_bounds:
                seqs_filtered[name] = seq
    if not seqs_filtered:
        raise ValueError("Too strict conditions")
    return seqs_filtered


def filter_length(seqs: dict, length_bounds=(0, 2**32)) -> dict:
    """
    Filters fastq dictionary by length.
    :param seqs: dict, fastq dictionary.
    :param length_bounds: float if one value is given (upper limit of filtration),
    tuple – otherwise (bounds of filtration). Default value is (0, 2**32).
    :return: dict, filtered fastq dictionary.
    Raises ValueError("Too strict conditions") if the return dictionary is empty.
    """
    seqs_filtered: dict = {}
    for name, seq in seqs.items():
        if isinstance(length_bounds, tuple):
            if length_bounds[0] <= len(seq[0]) <= length_bounds[1]:
                seqs_filtered[name] = seq
        else:
            if len(seq[0]) <= length_bounds:
                seqs_filtered[name] = seq
    if not seqs_filtered:
        raise ValueError("Too strict conditions")
    return seqs_filtered


def filter_quality(seqs: dict, quality_threshold=0) -> dict:
    """
    Filters fastq dictionary by quality.
    :param seqs: dict, fastq dictionary.
    :param quality_threshold: float, lower limit for filtration. Default value is 0.
    :return: dict, filtered fastq dictionary.
    Raises ValueError("Too strict conditions") if the return dictionary is empty.
    """
    seqs_filtered: dict = {}
    for name, seq in seqs.items():
        if compute_nucleotide_quality(seq[1]) >= quality_threshold:
            seqs_filtered[name] = seq
    if not seqs_filtered:
        raise ValueError("Too strict conditions")
    return seqs_filtered


def read_fastq_file(input_path: str) -> dict:
    """
    Reads a Fastq file and converts it into a dictionary sequentially.
    :param input_path: str, path to the input Fastq file.
    :return: dict, Fastq data in dictionary format.
    """
    fastq_data = {}

    with open(input_path, 'r') as file:
        header = False
        sequence = False
        comment = False
        quality = False

        for line in file:
            line = line.strip()
            if not header:
                header = line.strip('@')
            elif not sequence:
                sequence = line
            elif not comment:
                comment = line
            elif not quality:
                quality = line
                fastq_data[header] = [sequence, quality]
                header = False
                sequence = False
                comment = False
                quality = False

    return fastq_data


def write_fastq(fastq_data: dict, output_filename: str):
    """
    Writes fastq dictionary into the fastq file
    :param fastq_data: dict, fastq data
    :param output_filename: str, name of the output fastq file.
    """
    with open(output_filename, 'w') as fastq_file:
        for key, value in fastq_data.items():
            name = key
            sequence = value[0]
            quality = value[1]
            fastq_file.write(f'@{name}\n')
            fastq_file.write(f'{sequence}\n')
            fastq_file.write('+\n')
            fastq_file.write(f'{quality}\n')


def run_fastq_tools(input_path: str, output_filename=None,
                    gc_bounds=None, length_bounds=(0, 2**32), quality_threshold=0):
    """
    Filters fastq dictionary by GC content, length, and quality.
    :param seqs: dict, fastq dictionary.
    :param gc_bounds: float if one value is given (upper limit of filtration),
    tuple – otherwise (bounds of filtration). If no arguments are given, default value is None, returns input dictionary.
    :param quality_threshold: float, lower limit for filtration. Default value is 0.
    :param length_bounds: float if one value is given (upper limit of filtration),
    tuple – otherwise (bounds of filtration). Default value is (0, 2**32).
    :return: dict, filtered fastq dictionary.
    Raises ValueError("Too strict conditions") if the return dictionary in any of the functions is empty.
    """
    filtered_fastq = filter_quality(filter_gc_content(filter_length(seqs, length_bounds=length_bounds),
                                                      gc_bounds=gc_bounds), quality_threshold=quality_threshold)
    return filtered_fastq

def compute_gc_content(seq: str) -> float:
    """
    Computes GC-content of the input sequence.
    :param seq: nucleotide sequence
    :return: float, result of computation
    """
    return ((seq.count('G') + seq.count('C')) / len(seq)) * 100


def filter_gc_content(seqs: dict, gc_bounds=None) -> dict:
    """
    Filters fastq dictionary by GC-content.
    :param seqs: dict, fastq dictionary.
    :param gc_bounds: int or tuple
    :return: dict. Filtered fastq dictionary. If gc_bounds is not defined, returns input dictionary.
    """
    if gc_bounds is None:
        return seqs
    seqs_filtered: dict = dict()
    for name, seq in seqs.items():
        if type(gc_bounds) == tuple:
            if gc_bounds[0] <= compute_gc_content(seq[0]) <= gc_bounds[1]:
                seqs_filtered[name] = seq
        else:
            if compute_gc_content(seq[0]) <= gc_bounds:
                seqs_filtered[name] = seq
    if not seqs_filtered:
        raise ValueError("Too strict conditions")
    else:
        return seqs_filtered

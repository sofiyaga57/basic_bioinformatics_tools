def compute_gc_content(seq: str) -> float:
    """
    Computes GC-content of the input sequence.
    :param seq: nucleotide sequence
    :return: float, result of computation
    """
    return ((seq.count('G') + seq.count('C'))/len(seq)) * 100
import os


def read_fasta_file(input_fasta: str) -> dict:
    """
    Reads fasta file and saves it to dictionary
    :param input_fasta: str, fasta file
    :return dict, fasta dictionary
    """
    with open(os.path.abspath(input_fasta), mode='r') as fasta_file:
        fasta_data = {}
        current_name = None
        for line in fasta_file:
            line.strip()
            if line.startswith('>'):
                if current_name is not None:
                    fasta_data[current_name] = ''.join(seq)
                current_name = line.strip('>').strip()
                seq = []
            else:
                seq.append(line.strip())
        fasta_data[current_name] = ''.join(seq)
    return fasta_data

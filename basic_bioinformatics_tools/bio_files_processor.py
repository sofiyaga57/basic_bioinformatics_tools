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


def write_fasta_file(fasta_data: dict, output_fasta: str):
    """"
    Writes fasta dictionary into fasta file.
    :param fasta_data: dict, fasta dictionary
    :param output_fasta: str, name or path of output fasta file.
    Returns fasta file with written fasta dictionary
    """
    with open(output_fasta, mode='w') as file:
        for seq_id, seq in fasta_data.items():
            file.write('>' + seq_id + '\n')
            file.write(seq + '\n')


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta=None):
    """
    Converts multiline fasta file to oneline fasta file
    :param input_fasta: str, name of input fasta file with .fasta
    :param output_fasta: str, name of output fasta file without .fasta
    Returns file
    """
    input_fasta = os.path.abspath(input_fasta)
    if output_fasta is None:
        output_fasta = os.path.splitext(os.path.basename(input_fasta))[0]
    output_path = os.path.join(output_fasta + '.fasta')

    fasta_data = read_fasta_file(input_fasta)
    write_fasta_file(fasta_data=fasta_data, output_fasta=output_path)


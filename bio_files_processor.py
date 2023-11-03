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
        seq = []
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
    Returns file, names of sequences are seq1, seq2, etc.
    """
    input_fasta = os.path.abspath(input_fasta)
    if output_fasta is None:
        output_fasta = os.path.splitext(os.path.basename(input_fasta))[0]
    output_path = os.path.join(output_fasta + '.fasta')

    fasta_data = read_fasta_file(input_fasta)
    write_fasta_file(fasta_data=fasta_data, output_fasta=output_path)


def select_genes_from_gbk_to_list(input_gbk: str) -> list:
    """
    Selects gene's name and its translation into list
    :param input_gbk: str, name of input gbk file
    Returns list
    """

    gene_data = []
    in_cds = False
    in_translation = False
    current_gene = 'unknown'

    with open(input_gbk, mode='r') as gbk_file:
        for line in gbk_file:
            line = line.strip()
            if line.startswith("CDS"):
                in_cds = True
            if in_cds and line.startswith("/gene"):
                current_gene = line.split('gene="')[-1].split('"')[0]
            if in_cds and line.startswith("/translation"):
                in_translation = True
                current_translation = line.split('"/translation="')[-1].split('"')[1]
                continue
            if in_cds and in_translation:
                if line.endswith('"'):
                    current_translation += line.split(' ')[-1].split('"')[0]
                    in_translation = False
                    in_cds = False
                    gene_data.append(current_gene)
                    gene_data.append(current_translation)
                    current_gene = 'unknown'
                else:
                    current_translation += line

    return gene_data


def record_gene_translations(gene_data: list, genes: list, output_fasta: str, before=1, after=1):
    """
    Records translations of genes nearby the gene of interest into fasta file
    :param gene_data: list, list of gene's name and its translations
    :param genes: list, genes of interest
    :param before: int, number of genes to be recorded before the gene of interest
    :param after: int, number of genes to be recorded after the gene of interest
    :param output_fasta: str, name of output fasta file
    Returns fasta file
    """

    seq_num = 0  # index of sequence in fasta file

    with open(output_fasta, mode='w') as fasta_file:
        for interest_gene in genes:
            for k in range(len(gene_data)):
                if interest_gene in gene_data[k]:
                    index = k
                    left_index = index - 2 * before
                    right_index = index + 2 * after  # define ranges for genes to write into fasta file
                    for i in range(left_index + 1, index, 2):
                        seq_num += 1
                        fasta_file.write(f">seq{seq_num}\n")
                        fasta_file.write(f"{gene_data[i]}\n")
                    for j in range(right_index + 1, index + 1, -2):
                        seq_num += 1
                        fasta_file.write(f">seq{seq_num}\n")
                        fasta_file.write(f"{gene_data[j]}\n")


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: list, n_before=1, n_after=1, output_fasta=None):
    """
    Selects genes surrounding genes of interest.
    :param input_gbk: file in a gbk format
    :param genes: list, genes of interest
    :param n_before: int, number of genes to be recorded before the gene of interest
    :param n_after: int, number of genes to be recorded after the gene of interest
    :param output_fasta: str, name of output fasta file without .fasta
    Returns fasta file
    """
    input_gbk = os.path.abspath(input_gbk)
    if output_fasta is None:
        output_fasta = os.path.splitext(os.path.basename(input_gbk))[0]
    output_path = os.path.join(output_fasta + '.fasta')

    gene_data = select_genes_from_gbk_to_list(input_gbk=input_gbk)
    record_gene_translations(gene_data=gene_data, genes=genes, before=n_before,
                             after=n_after, output_fasta=output_path)


def shift_sequence(seq: str, shift: int) -> str:
    """
    Records translations of genes nearby the gene of interest into fasta file
    :param seq: str, sequence to shift
    :param  shift: int, how many nucleotides to shift
    Returns str, shifted sequence
    """
    seq_list = list(seq)
    shift = shift % len(seq)
    if shift > 0:
        for i in range(shift):
            seq_list.insert(0, seq_list.pop())
    else:
        shift = abs(shift)
        for i in range(shift):
            seq_list.insert(len(seq_list)-1, seq_list.pop(0))
    return ''.join(seq_list)


def change_fasta_start_pos(input_fasta: str, shift: int, output_fasta=None):
    """
    Records translations of genes nearby the gene of interest into fasta file
    :param input_fasta: str, input fasta file
    :param shift: int, how many nucleotides to shift
    :param output_fasta; str, name of output_file, without .fasta, optional.
    If not mentioned, input_fasta name is used.
    Return fasta file with shifted sequence.
    """

    shifted_fasta_data = {}

    input_fasta = os.path.abspath(input_fasta)
    if output_fasta is None:
        output_fasta = os.path.splitext(os.path.basename(input_fasta))[0]
    output_path = os.path.join(f'{output_fasta}_shifted.fasta')

    fasta_data = read_fasta_file(input_fasta=input_fasta)
    sequence, name = list(fasta_data.values())[0], list(fasta_data.keys())[0]
    shifted_seq = shift_sequence(seq=sequence, shift=shift)
    shifted_fasta_data[name] = shifted_seq
    write_fasta_file(fasta_data=shifted_fasta_data, output_fasta=output_path)

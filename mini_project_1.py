from Bio import SeqIO

def sequence_extraction(fasta_file, contig_name, start, end):
    for record in SeqIO.parse(fasta_file_name, "fasta"):
        if record.id == num_sequence:
            sequence = str(record.seq)
            start_index = start - 1
            end_index = end
            extracted_sequence = sequence[start_index:end_index]
            return extracted_sequence
    return None


if __name__ == "__main__":
    # Input FASTA file, contig name, and coordinates
    fasta_file_name = input("Input FASTA:")
    #fasta_file_name = 'example.fasta'
    num_sequence = input("Contig Name:")
    start = int(input("Start Coordinate:"))
    end = int(input("End Coordinate:"))
    result = sequence_extraction(fasta_file_name, num_sequence, start, end)
    print (result)









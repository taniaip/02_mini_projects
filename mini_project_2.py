import argparse
import gffutils
from Bio import SeqIO


# Set up argument parser
parser = argparse.ArgumentParser()

parser.add_argument(
        "-f", "--fasta",
        dest='fasta_file',
        type=str,
        required=True,
        help="Input the fasta file")

parser.add_argument(
        "-c", "--contig",
        dest='contig_name',
        type=str,
        required=True,
        help="Input the contig name ")

parser.add_argument(
        "-s", "--start",
        dest='start',
        type=int,
        required=True,
        help="Input the start")

parser.add_argument(
        "-e", "--end",
        dest='end',
        type=int,
        required=True,
        help="Input the end")


parser.add_argument(
        "-g", "--gff3",
        dest='gff3_file',
        type=str,
        required=True,
        help="Input GFF3 genome file"


# Parse the arguments
args = parser.parse_args()

db = gffutils.create_db('stop_motifs_updtaed.gff3', ':memory:')


def sequence_extraction(fasta_file, contig_name, start, end):
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == contig_name:
            sequence = str(record.seq)
            start_index = start - 1  # Python uses 0-based indexing
            end_index = end  # 'end' is inclusive in Python slicing
            extracted_sequence = sequence[start_index:end_index]
            return extracted_sequence
    return None

def main(args):
    # Call the function and print the result
    result = sequence_extraction(args.fasta_file, args.contig_name, args.start, args.end)
    print(result)


if __name__ == '__main__':
    main(args)

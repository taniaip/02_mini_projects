import argparse
import gffutils
from pprint import pprint
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
        "-g", "--gff3",
        dest='gff3_file',
        type=str,
        required=True,
        help="Input GFF3 genome file")


# Parse the arguments
args = parser.parse_args()


def list_extraction_motif(gff3_file):
    db = gffutils.create_db(args.gff3_file, dbfn=':memory:')
    return [ (motif.seqid, motif.start, motif.end) for motif in db.features_of_type('motif') ]

def sequence_extraction(fasta_file, location_list):
    seq_record= SeqIO.parse(args.fasta_file, "fasta"):
    for record in seq_record:
        extracted_sequences = []
        for motif in location_list:
            contig_name, start, end = motif
            if contig_name == record.id:
                sequence = str(record.seq)
                extracted_sequences.append( sequence[start-1:end] )
            
    return extracted_sequences


def main(args):
    result_location = list_extraction_motif(args.gff3_file)
    pprint(sequence_extraction(args.fasta_file, result_location)):


if __name__ == '__main__':
    main(args)    
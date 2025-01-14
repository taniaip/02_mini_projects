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
    fasta_sequences = {record.id : record.seq for record in SeqIO.parse(fasta_file, "fasta")}    
    extracted_sequences = []
    for motif in location_list:
        contig_name, start, end = motif
        for record in fasta_sequences:
            if contig_name == record:
                sequence = str(record.seq)
                extracted_sequences.append( sequence[start-1:end])
            
    return extracted_sequences


def main(args):
    result_location = list_extraction_motif(args.gff3_file)
    result_sequence = sequence_extraction(args.fasta_file, result_location)
    print(result_sequence)

if __name__ == '__main__':
    main(args)    











# oct 23rd project
# import argparse
# import gffutils
# from pprint import pprint
# from Bio import SeqIO
# from Bio.Seq import Seq


# # Set up argument parser
# parser = argparse.ArgumentParser()

# parser.add_argument(
#         "-f", "--fasta",
#         dest='fasta_file',
#         type=str,
#         required=True,
#         help="Input the fasta file")

# parser.add_argument(
#         "-g", "--gff3",
#         dest='gff3_file',
#         type=str,
#         required=True,
#         help="Input GFF3 genome file")

# args = parser.parse_args()


# def list_polyA_location(gff3_file):
#     db = gffutils.create_db(args.gff3_file, dbfn=':memory:')
#     return [ (motif.seqid, motif.start, motif.strand) for motif in db.features_of_type('polyAsite') ]

# def sequence_extraction(fasta_file, location_list):
#     fasta_sequences = {record.id: record.seq for record in SeqIO.parse(fasta_file, "fasta")}    
#     extracted_sequences = []
#     for motif in location_list:
#         contig_name, polyasite, strandType= polyAsite
#         for record in fasta_sequences:
#             if contig_name == record:
#                 sequence = str(record.seq)
#                 if strandType == '+':
#                     extracted_sequences.append( sequence[polyasite-25:polyasite+25])
#                 else:
#                     extracted_sequences.append( str(record.seq[polyasite-25:polyasite+25].reverse_complement()))

            
#     return extracted_sequences

# def main(args):
#     result_location = list_polyA_location(args.gff3_file)
#     result_sequence = sequence_extraction(args.fasta_file, result_location)
#     print(result_sequence)

# if __name__ == '__main__':
#     main(args)        
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

    # extracted_motifs = []
    # for motif in db.all_features():
    # #for motif in db.features_of_type('motif'):
    #     location = (motif.seqid, motif.start, motif.end)
    #     extracted_motifs.append(location)
    # return extracted_motifs

    # extracted_motifs = []
    # for motif in db.features_of_type('motif'):
    #     extracted_motifs.append((motif.seqid, motif.start, motif.end))

    # extracted_motifs = [ (motif.seqid, motif.start, motif.end) for motif in db.features_of_type('motif') ]
    # return extracted_motifs



# def sequence_extraction(fasta_file, location_list):
    # extracted_sequence = []
    # for location in location_list:
    #     contig_name, start, end = location
    #     for record in SeqIO.parse(fasta_file, "fasta"):
    #         if contig_name == record.id:
    #             sequence = str(record.seq)
    #             extracted_sequence.append(sequence[start-1:end])
    # return extracted_sequence
    # return "nothing is found"   

def sequence_extraction(record, location_list):
    extracted_sequences = []
    for motif in location_list:
        contig_name, start, end = motif
        if contig_name == record.id:
            sequence = str(record.seq)
            extracted_sequences.append( sequence[start-1:end] )
            
    return extracted_sequences


def main(args):
    result_location = list_extraction_motif(args.gff3_file)
    #pprint(result_location)
    #print(len(result_location))
    print("hi")

    for record in SeqIO.parse(args.fasta_file, "fasta"):
        result_sequences = sequence_extraction(record, result_location)
        pprint(result_sequences)
        # for contig_name, start, end in result_location:
        #     if contig_name == record.id:
        #         sequence = str(record.seq)
        #         extracted_sequence = sequence[start-1:end]
        #         print(extracted_sequence)


    #result_sequence = sequence_extraction(args.fasta_file, result_location)    
    print("hi")
    #print(headresult_sequence)
    print(len(container))


if __name__ == '__main__':
    main(args)    
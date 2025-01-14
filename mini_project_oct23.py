import argparse
import gffutils
from pprint import pprint
from Bio import SeqIO
from Bio.Seq import Seq


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

args = parser.parse_args()


def list_polyA_location(gff3_file):
    db = gffutils.create_db(args.gff3_file, dbfn=':memory:', merge_strategy='create_unique')
    return [ (polyA.seqid, polyA.start, polyA.strand) for polyA in db.features_of_type('polyA') ]


def sequence_extraction(fasta_file, location_list):
    # fasta_sequences = {record.id : record for record in SeqIO.parse(fasta_file, "fasta")}
    
    # collect all sequence records in a list
    fasta_sequences = [record for record in SeqIO.parse(fasta_file, "fasta")]

    # capture "all polyA site region" sequences
    extracted_sequences = []
    # compare all polyA sites with all seqrecords
    for polyAsites_info in location_list:
        contig_name, polyAsite, strandType = polyAsites_info
        for record in fasta_sequences:

            # only extract sequences if polyA site contig
            # matches seq record contig
            if contig_name == record.id:
                sequence = str(record.seq)

                polyA_region_seq = [
                    sequence[polyAsite-26:polyAsite-1].upper(),
                    sequence[polyAsite-1].upper(),
                    sequence[polyAsite:polyAsite+25].upper()
                    # sequence[polyAsite-25:polyAsite+25].upper()
                ]

                if strandType == '+':
                    extracted_sequences.append(polyA_region_seq)
                # elif strandType == '-':
                #     polyA_region_seq_rc = [  str( Seq(seq).reverse_complement() ) for seq in polyA_region_seq.reverse() ]
                #     # (Seq(polyA_region_seq)).reverse_complement()
                #     # print(contig_name, polyAsite, strandType, polyA_region_seq_rc)
                #     # extracted_sequences.append( str(polyA_region_seq_rc) )
                #     extracted_sequences.append( polyA_region_seq_rc )
                break    
    return extracted_sequences


def main(args):
    location_list = list_polyA_location(args.gff3_file)
    # pprint(location_list)
    result_sequence = sequence_extraction(args.fasta_file, location_list)
    pprint(result_sequence)


if __name__ == '__main__':
    main(args)        



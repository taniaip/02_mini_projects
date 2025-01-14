import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import weblogo as wl
import io

# Argument parser setup
parser = argparse.ArgumentParser(description="Generate sequence logos from GFF3 and FASTA files.")
parser.add_argument(
    "-g", "--gff3",
    dest="gff3_file",
    type=str,
    required=True,
    help="Input GFF3 file with polyA site information."
)
parser.add_argument(
    "-f", "--fasta",
    dest="fasta_file",
    type=str,
    required=True,
    help="Input FASTA file containing genomic sequences."
)
parser.add_argument(
    "-o", "--out",
    dest="out_prefix",
    type=str,
    required=True,
    help="Output prefix for the sequence logo (PNG format)."
)


def extract_polyA_sequences(gff3_file, fasta_file, flank=24):
    """
    Extract sequences around polyA sites specified in the GFF3 file from the FASTA file.
    :param gff3_file: GFF3 file with polyA site information.
    :param fasta_file: FASTA file containing genomic sequences.
    :param flank: Number of bases to extract upstream and downstream of the polyA site.
    :return: List of sequences around polyA sites.
    """
    # Load sequences from FASTA
    fasta_sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

    # Extract sequences from GFF3
    sequences = []
    with open(gff3_file, "r") as gff3:
        for line in gff3:
            if line.startswith("#"):
                continue
            columns = line.strip().split("\t")
            if len(columns) < 9 or columns[2] != "polyA":  # Ensure polyA feature type
                continue

            seqid = columns[0]  # Contig/Sequence ID
            start = int(columns[3]) - 1  # 0-based start position
            strand = columns[6]  # Strand

            if seqid not in fasta_sequences:
                continue  # Skip if contig is not in the FASTA file

            # Extract the neighborhood sequence
            contig_sequence = fasta_sequences[seqid]
            polyA_region_seq = contig_sequence[max(0, start - flank):start + flank + 1].upper()

            # Reverse complement if the strand is negative
            if strand == "-":
                polyA_region_seq = str(Seq(polyA_region_seq).reverse_complement())

            sequences.append(polyA_region_seq)

    return sequences


def generate_logo(sequences, output_file):
    """
    Generate a sequence logo from a list of sequences and save it as an image.
    :param sequences: List of sequences to generate the logo.
    :param output_file: Output file name for the logo image (PNG format).
    """
    # Convert sequences to a single newline-separated string
    seq_data = "\n".join(sequences)
    seqs = wl.read_seq_data(io.StringIO(seq_data))

    # Create logo data
    logo_data = wl.LogoData.from_seqs(seqs)

    # Set up options for the sequence logo
    logo_options = wl.LogoOptions()
    logo_options.title = "Sequence Logo"
    logo_options.color_scheme = wl.std_color_schemes['classic']
    logo_options.logo_width = 60 # Set the width based on sequence length
    logo_options.stack_width = 12  # Adjust the width of each stack

    # Format and generate the sequence logo
    logo_format = wl.LogoFormat(logo_data, logo_options)
    png = wl.png_print_formatter(logo_data, logo_format)

    # Save the PNG file
    with open(output_file, "wb") as out:
        out.write(png)


def main(args):
    # Extract sequences from the GFF3 and FASTA files
    sequences = extract_polyA_sequences(args.gff3_file, args.fasta_file)

    # Generate and save the sequence logo
    generate_logo(sequences, f"{args.out_prefix}.png")


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)










# #!/usr/bin/env python

# # standard modules
# import argparse
# import itertools
# from pprint import pprint

# # modules you'll need to install
# import gffutils
# from Bio import SeqIO
# from Bio.Seq import Seq


# # Set up argument parser
# parser = argparse.ArgumentParser()

# parser.add_argument(
#         "-f", "--fasta",
#         dest='fasta_file',
#         type=str,
#         required=True,
#         help="Input the FASTA file")

# parser.add_argument(
#         "-g", "--gff3",
#         dest='gff3_file',
#         type=str,
#         required=True,
#         help="Input GFF3 file with polyA sites")

# args = parser.parse_args()


# def main(args):

#     # load GFF3 file
#     db = gffutils.create_db(args.gff3_file, dbfn=':memory:', merge_strategy='create_unique')

#     # get polyA site locations (contig, start, end)
#     polyA_locations = [ 
#         (polyA.id, polyA.seqid, polyA.start, polyA.strand)
#         for polyA in db.features_of_type('polyA')
#     ]

#     # load contigs in an iterator of SeqRecord objects
#     fasta_records = SeqIO.parse(args.fasta_file, "fasta")

#     # compare all polyA sites with all SeqRecords
#     for location, record in itertools.product(polyA_locations, fasta_records):

#         # skip combinations where polyA contig and FASTA contig do not match
#         if location[1] != record.id: continue

#         # get the sequence of the polyA site neighborhood
#         polyA_region_seq = extract_polyA_neighborhood_seq(location, record)

#         # break 'neighborhood' down in sections
#         upstream   = polyA_region_seq[:24]
#         # upstream   = polyA_region_seq[15:24]
#         polyA_nt   = polyA_region_seq[24] # the polyA site coordinate
#         buffer     = polyA_region_seq[25:29]
#         motif      = polyA_region_seq[29:37] # the TGTTTGTT motif
#         downstream = polyA_region_seq[37:]

#         # print()
#         # print(f'{location = }\t{polyA_region_seq}')

#         # print(f'{location = }\t{upstream}--{polyA_nt}-{buffer}-{motif}--{downstream}')

#         print(f">{location[1]}_{location[2]}")
#         print(f'{upstream}{polyA_nt}{buffer}{motif}{downstream}')



# def extract_polyA_neighborhood_seq(location, record):

#     # unpacking location and getting contig seq
#     polyA_id, polyA_contig, polyA_site, polyA_strand = location
#     contig_sequence = str(record.seq)

#     # convert polyA_site (1-indexed) to 0-indexed
#     polyA_site = polyA_site - 1

#     # extract polyA neighborhood seq,
#     # such that the polyA site is EXACTLY in the middle!
#     ## this ensures that after reverse complementing the sequence for - strand polyA sites,
#     ## the polyA site is still in the same position within the neighborhood as for + strand polyA sites
#     polyA_region_seq = contig_sequence[polyA_site-24:polyA_site+25].upper()
#     if polyA_strand == '-':
#         polyA_region_seq = Seq(polyA_region_seq).reverse_complement()

#     return polyA_region_seq




# if __name__ == '__main__':
#     main(args)





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
    logo_options.logo_width = 50  # set the width to 50 characters  #but it didn't work!!:)
    logo_options.stack_width = 12  

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

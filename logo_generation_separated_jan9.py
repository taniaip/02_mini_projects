import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import weblogo as wl
import regex 
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
    # Initialize category lists
    sequences_3 = []
    sequences_4 = []
    sequences_5 = []
    sequences_else = []
    sequences = []

    # Precompile your three patterns.
    # Explanation: 
    #   .{24}   --> skip first 24 positions
    #   ([ACGT]{3}) --> exactly 3 bases
    #   ([ACGT]{8}){e<=2} --> 8-base substring with up to 2 mismatches to "TGTTTGTT"

    # We'll do a simpler approach: we match the explicit "TGTTTGTT" with {s<=2} for substitutions only, ignoring insertions/deletions.
    pattern_3 = regex.compile(rf"^(.{{24}})([ACGT]{{3}})(TGTTTGTT){{s<=2}}", flags=regex.IGNORECASE)
    pattern_4 = regex.compile(rf"^(.{{24}})([ACGT]{{4}})(TGTTTGTT){{s<=2}}", flags=regex.IGNORECASE)
    pattern_5 = regex.compile(rf"^(.{{24}})([ACGT]{{5}})(TGTTTGTT){{s<=2}}", flags=regex.IGNORECASE)

    # Load FASTA
    fasta_seqs = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

    # Parse GFF3
    with open(gff3_file, "r") as gf:
        for line in gf:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "polyA":
                continue

            seqid = parts[0]
            start = int(parts[3]) - 1  # 0-based
            strand = parts[6]

            if seqid not in fasta_seqs:
                continue

            full_seq = fasta_seqs[seqid]
            # Extract ± flank
            region_seq = full_seq[max(0, start - flank) : start + flank + 1].upper()

            # Reverse complement if negative
            if strand == "-":
                region_seq = str(Seq(region_seq).reverse_complement())
            
            sequences.append(region_seq)

            # Attempt to match with each pattern
            matched_3 = pattern_3.match(region_seq)
            matched_4 = pattern_4.match(region_seq)
            matched_5 = pattern_5.match(region_seq)

            if matched_3:
                sequences_3.append(region_seq)
            elif matched_4:
                sequences_4.append(region_seq)
            elif matched_5:
                sequences_5.append(region_seq)
            else:
                sequences_else.append(region_seq)

    print("Group with 3-gap and TGTTTGTT (<=2 subs):", len(sequences_3))
    print("Group with 4-gap:", len(sequences_4))
    print("Group with 5-gap:", len(sequences_5))
    print("Others (no match):", len(sequences_else))

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

#!/usr/bin/env python

###For has_stop category:
##(base) tania@tanias-MacBook-Air 02_mini_projects % python logo_generation_separated_jan9.py -g1 march4_has_stop_polyA.gff3 
# -f ST2_sorted_masked.fasta --outAll has_stop_all --out3 has_stop_gap3 --out4a has_stop_gap4a --out4b has_stop_gap4b 
# --out4c has_stop_gap4c --out5 has_stop_gap5 --outElse has_stop_else -g2 else_polyAs_has_stop_logo.gff3

###Results for has_stop category:
# Total sequences: 5163
# 3-gap matches: 52
# 4-gap matches, just T: 1828
# 4-gap matches, TG: 372
# 4-gap matches, TA: 215
# 5-gap matches: 175
# else (no match): 2521


import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import regex
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
import io

def parse_args():
    parser = argparse.ArgumentParser(description="Generate sequence logos (PNG) and write 'else' GFF3 lines.")
    parser.add_argument("-g1", "--gff3_in", dest="gff3_in", type=str, required=True,
                        help="Input GFF3 file with polyA site info.")
    parser.add_argument("-g2", "--gff3_else", dest="gff3_else", type=str, required=True,
                        help="Output GFF3 for 'else' category.")
    parser.add_argument("-f", "--fasta", dest="fasta_file", type=str, required=True,
                        help="Input FASTA file of genomic sequences.")
    parser.add_argument("--outAll",  dest="out_all",  required=True, help="Logo for all sequences.")
    parser.add_argument("--out3",    dest="out_3",    required=True, help="Logo for 3-gap sequences.")
    parser.add_argument("--out4a",   dest="out_4a",   required=True, help="Logo for 4-gap (just T).")
    parser.add_argument("--out4b",   dest="out_4b",   required=True, help="Logo for 4-gap (TG).")
    parser.add_argument("--out4c",   dest="out_4c",   required=True, help="Logo for 4-gap (TA).")
    parser.add_argument("--out5",    dest="out_5",    required=True, help="Logo for 5-gap sequences.")
    parser.add_argument("--outElse", dest="out_else", required=True, help="Logo for 'else' sequences.")
    return parser.parse_args()

def extract_polyA_sequences(gff3_in, fasta_file, flank=24):
    sequences_all  = []
    sequences_3    = []
    sequences_4a   = []
    sequences_4b   = []
    sequences_4c   = []
    sequences_5    = []
    sequences_else = []
    else_gff_lines = []

    pattern_3  = regex.compile(r"^(.{24})(T)([ACGT]{3})(TGTTTGTT){s<=2}")
    pattern_4a = regex.compile(r"^(.{24})(T)([ACGT]{4})(TGTTTGTT){s<=2}")
    pattern_4b = regex.compile(r"^(.{23})(TG)([ACGT]{4})(TGTTTGTT){s<=2}")
    pattern_4c = regex.compile(r"^(.{23})(TA)([ACGT]{4})(TGTTTGTT){s<=2}")
    pattern_5  = regex.compile(r"^(.{24})(T)([ACGT]{5})(TGTTTGTT){s<=2}")

    fasta_seqs = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

    with open(gff3_in, "r") as gf:
        for line in gf:
            line = line.strip()
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9 or parts[2] != "polyA":
                continue

            seqid = parts[0]
            start = int(parts[3]) - 1
            strand = parts[6]
            if seqid not in fasta_seqs:
                continue

            full_seq = fasta_seqs[seqid]
            region_seq = full_seq[max(0, start - flank): start + flank + 1].upper()
            if strand == "-":
                region_seq = str(Seq(region_seq).reverse_complement())

            sequences_all.append(region_seq)
            matched_3  = pattern_3.match(region_seq)
            matched_4a = pattern_4a.match(region_seq)
            matched_4b = pattern_4b.match(region_seq)
            matched_4c = pattern_4c.match(region_seq)
            matched_5  = pattern_5.match(region_seq)

            if matched_4a:
                sequences_4a.append(region_seq)
            elif matched_4b:
                sequences_4b.append(region_seq)
            elif matched_4c:
                sequences_4c.append(region_seq)
            elif matched_3:
                sequences_3.append(region_seq)
            elif matched_5:
                sequences_5.append(region_seq)
            else:
                sequences_else.append(region_seq)
                else_gff_lines.append(line)

    print(f"Total sequences: {len(sequences_all)}")
    print(f"3-gap matches: {len(sequences_3)}")
    print(f"4-gap matches, just T: {len(sequences_4a)}")
    print(f"4-gap matches, TG: {len(sequences_4b)}")
    print(f"4-gap matches, TA: {len(sequences_4c)}")
    print(f"5-gap matches: {len(sequences_5)}")
    print(f"else (no match): {len(sequences_else)}")
    return sequences_all, sequences_3, sequences_4a, sequences_4b, sequences_4c, sequences_5, sequences_else, else_gff_lines

def save_logo_logomaker(sequences, output_file):
    if not sequences:
        print(f"No sequences found for {output_file}, skipping logo.")
        return
    # Convert the alignment to an information content matrix (bits)
    info_df = logomaker.alignment_to_matrix(sequences, to_type='information')

    # Set up the logo figure
    plt.figure(figsize=(len(info_df)/3, 3))
    logo = logomaker.Logo(info_df, color_scheme='classic')
    logo.style_spines(visible=False)
    logo.style_spines(spines=['left', 'bottom'], visible=True)
    logo.ax.set_ylabel('bits')
    logo.ax.set_ylim([0, 2.0])  # Y-axis from 0 to 2 bits
    plt.title(f"Sequence Logo: {output_file}")
    plt.tight_layout()
    # Save as PNG
    plt.savefig(output_file + ".png", dpi=200)
    plt.close()

def main():
    args = parse_args()
    (seq_all, seq_3,seq_4a,
     seq_4b, seq_4c, seq_5, seq_else, else_gff_lines) = extract_polyA_sequences(
        gff3_in=args.gff3_in,
        fasta_file=args.fasta_file,
        flank=24
    )

    # Save logos as PNG (will be named e.g., outAll.png, out3.png, ...)
    save_logo_logomaker(seq_all, args.out_all)
    save_logo_logomaker(seq_3, args.out_3)
    save_logo_logomaker(seq_4a, args.out_4a)
    save_logo_logomaker(seq_4b, args.out_4b)
    save_logo_logomaker(seq_4c, args.out_4c)
    save_logo_logomaker(seq_5, args.out_5)
    save_logo_logomaker(seq_else, args.out_else)

    # Write else GFF3
    with open(args.gff3_else, "w") as f:
        for gff_line in else_gff_lines:
            f.write(gff_line + "\n")

if __name__ == "__main__":
    main()


#### This script worked on the lab's computer but not on my laptop, but the new version works on my laptop:
# # $ python logo_generation_separated_jan9.py -g1 march4_non_stop_polyA.gff3 -f ST2_sorted_masked.fasta --outAll non_stop_all --out3 non_stop_gap3 --out4a non_stop_gap4a --out4b non_stop_gap4b --out4c non_stop_gap4c --out5 non_stop_gap5 --outElse non_stop_else -g2 else_polyAs_logo.gff3
# # Total sequences: 968
# # 3-gap matches: 3
# # 4-gap matches, just T: 673
# # 4-gap matches, TG: 86
# # 4-gap matches, TA: 88
# # 5-gap matches: 10
# # else (no match): 108

# #!/usr/bin/env python
# import argparse
# from Bio import SeqIO
# from Bio.Seq import Seq
# import weblogo as wl
# import regex
# import io

# def parse_args():
#     parser = argparse.ArgumentParser(
#         description="Generate separate sequence logos and save 'else' category lines to a separate GFF3."
#     )
#     # GFF3 input (main polyA sites)
#     parser.add_argument(
#         "-g1", "--gff3_in",
#         dest="gff3_in",
#         type=str,
#         required=True,
#         help="Input GFF3 file with polyA site information."
#     )
#     # GFF3 output (for 'else' category lines)
#     parser.add_argument(
#         "-g2", "--gff3_else",
#         dest="gff3_else",
#         type=str,
#         required=True,
#         help="Output GFF3 file to store lines that don't match any pattern."
#     )
#     # FASTA input
#     parser.add_argument(
#         "-f", "--fasta",
#         dest="fasta_file",
#         type=str,
#         required=True,
#         help="Input FASTA file containing genomic sequences."
#     )
#     # Output logos
#     parser.add_argument("--outAll",  dest="out_all",  required=True, help="Logo for all sequences.")
#     parser.add_argument("--out3",    dest="out_3",    required=True, help="Logo for 3-gap sequences.")
#     parser.add_argument("--out4a",   dest="out_4a",   required=True, help="Logo for 4-gap (just T).")
#     parser.add_argument("--out4b",   dest="out_4b",   required=True, help="Logo for 4-gap (TG).")
#     parser.add_argument("--out4c",   dest="out_4c",   required=True, help="Logo for 4-gap (TA).")
#     parser.add_argument("--out5",    dest="out_5",    required=True, help="Logo for 5-gap sequences.")
#     parser.add_argument("--outElse", dest="out_else", required=True, help="Logo for 'else' category sequences.")

#     return parser.parse_args()

# def extract_polyA_sequences(gff3_in, fasta_file, flank=24):
#     """
#     Extract sequences around polyA sites, classify them by various fuzzy patterns,
#     and store lines that match "else" category in a separate list for writing later.
#     """
#     sequences_all  = []
#     sequences_3    = []
#     sequences_4a   = []
#     sequences_4b   = []
#     sequences_4c   = []
#     sequences_5    = []
#     sequences_else = []
#     else_gff_lines = []

#     # Patterns for 3-gap, 4-gap, 5-gap
#     pattern_3  = regex.compile(r"^(.{24})(T)([ACGT]{3})(TGTTTGTT){s<=2}")
#     pattern_4a = regex.compile(r"^(.{24})(T)([ACGT]{4})(TGTTTGTT){s<=2}")
#     pattern_4b = regex.compile(r"^(.{23})(TG)([ACGT]{4})(TGTTTGTT){s<=2}")
#     pattern_4c = regex.compile(r"^(.{23})(TA)([ACGT]{4})(TGTTTGTT){s<=2}")
#     pattern_5  = regex.compile(r"^(.{24})(T)([ACGT]{5})(TGTTTGTT){s<=2}")

#     # Load FASTA
#     fasta_seqs = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

#     with open(gff3_in, "r") as gf:
#         for line in gf:
#             line = line.strip()
#             if line.startswith("#"):
#                 continue
#             parts = line.split("\t")
#             if len(parts) < 9 or parts[2] != "polyA":
#                 continue

#             seqid = parts[0]
#             start = int(parts[3]) - 1
#             strand = parts[6]

#             if seqid not in fasta_seqs:
#                 continue

#             full_seq = fasta_seqs[seqid]
#             region_seq = full_seq[max(0, start - flank) : start + flank + 1].upper()

#             # Reverse complement if negative
#             if strand == "-":
#                 region_seq = str(Seq(region_seq).reverse_complement())

#             sequences_all.append(region_seq)

#             # Check patterns
#             matched_3  = pattern_3.match(region_seq)
#             matched_4a = pattern_4a.match(region_seq)
#             matched_4b = pattern_4b.match(region_seq)
#             matched_4c = pattern_4c.match(region_seq)
#             matched_5  = pattern_5.match(region_seq)

#             if matched_4a:
#                 sequences_4a.append(region_seq)
#             elif matched_4b:
#                 sequences_4b.append(region_seq)
#             elif matched_4c:
#                 sequences_4c.append(region_seq)
#             elif matched_3:
#                 sequences_3.append(region_seq)
#             elif matched_5:
#                 sequences_5.append(region_seq)
#             else:
#                 sequences_else.append(region_seq)
#                 else_gff_lines.append(line)

#     # Print summary
#     print(f"Total sequences: {len(sequences_all)}")
#     print(f"3-gap matches: {len(sequences_3)}")
#     print(f"4-gap matches, just T: {len(sequences_4a)}")
#     print(f"4-gap matches, TG: {len(sequences_4b)}")
#     print(f"4-gap matches, TA: {len(sequences_4c)}")
#     print(f"5-gap matches: {len(sequences_5)}")
#     print(f"else (no match): {len(sequences_else)}")

#     return sequences_all, sequences_3, sequences_4a, sequences_4b, sequences_4c, sequences_5, sequences_else, else_gff_lines

# def generate_logo(sequences, output_file):
#     """
#     Generate a sequence logo from a list of sequences and save it as an image.
#     """
#     if not sequences:
#         print(f"No sequences found for {output_file}, skipping logo.")
#         return

#     seq_data = "\n".join(sequences)
#     seqs = wl.read_seq_data(io.StringIO(seq_data))

#     logo_data = wl.LogoData.from_seqs(seqs)
#     logo_options = wl.LogoOptions()
#     logo_options.title = f"Sequence Logo: {output_file}"
#     logo_options.color_scheme = wl.std_color_schemes['classic']

#     logo_format = wl.LogoFormat(logo_data, logo_options)
#     png = wl.png_print_formatter(logo_data, logo_format)

#     with open(output_file, "wb") as out:
#         out.write(png)

# def main():
#     args = parse_args()

#     # 1. Extract & classify sequences
#     (seq_all, seq_3,seq_4a,
#      seq_4b, seq_4c, seq_5, seq_else, else_gff_lines) = extract_polyA_sequences(
#         gff3_in=args.gff3_in,
#         fasta_file=args.fasta_file,
#         flank=24
#     )

#     # 2. Generate logos for each group
#     generate_logo(seq_all, f"{args.out_all}.png")
#     generate_logo(seq_3, f"{args.out_3}.png")
#     generate_logo(seq_4a, f"{args.out_4a}.png")
#     generate_logo(seq_4b, f"{args.out_4b}.png")
#     generate_logo(seq_4c, f"{args.out_4c}.png")
#     generate_logo(seq_5, f"{args.out_5}.png")
#     generate_logo(seq_else, f"{args.out_else}.png")

#     # 3. Write "else" category lines to an output GFF3
#     with open(args.gff3_else, "w") as f:
#         for gff_line in else_gff_lines:
#             f.write(gff_line + "\n")

# if __name__ == "__main__":
#     main()
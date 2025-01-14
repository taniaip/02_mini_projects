# # #!/usr/bin/env python

# # # standard modules
# # import argparse
# # import itertools
# # from pprint import pprint

# # # modules you'll need to install
# # import gffutils
# # from Bio import SeqIO
# # from Bio.Seq import Seq


# # # Set up argument parser
# # parser = argparse.ArgumentParser()

# # parser.add_argument(
# #         "-f", "--fasta",
# #         dest='fasta_file',
# #         type=str,
# #         required=True,
# #         help="Input the FASTA file")

# # parser.add_argument(
# #         "-g", "--gff3",
# #         dest='gff3_file',
# #         type=str,
# #         required=True,
# #         help="Input GFF3 file with polyA sites")

# # args = parser.parse_args()


# # def main(args):

# #     # load GFF3 file
# #     db = gffutils.create_db(args.gff3_file, dbfn=':memory:', merge_strategy='create_unique')

# #     # get polyA site locations (contig, start, end)
# #     polyA_locations = [ 
# #         (polyA.id, polyA.seqid, polyA.start, polyA.strand)
# #         for polyA in db.features_of_type('polyA')
# #     ]

# #     # load contigs in an iterator of SeqRecord objects
# #     fasta_records = SeqIO.parse(args.fasta_file, "fasta")

# #     # compare all polyA sites with all SeqRecords
# #     for location, record in itertools.product(polyA_locations, fasta_records):

# #         # skip combinations where polyA contig and FASTA contig do not match
# #         if location[1] != record.id: continue

# #         # get the sequence of the polyA site neighborhood
# #         polyA_region_seq = extract_polyA_neighborhood_seq(location, record)

# #         # break 'neighborhood' down in sections
# #         upstream   = polyA_region_seq[:24]
# #         # upstream   = polyA_region_seq[15:24]
# #         polyA_nt   = polyA_region_seq[24] # the polyA site coordinate
# #         buffer     = polyA_region_seq[25:29]
# #         motif      = polyA_region_seq[29:37] # the TGTTTGTT motif
# #         downstream = polyA_region_seq[37:]

# #         # print()
# #         # print(f'{location = }\t{polyA_region_seq}')

# #         # print(f'{location = }\t{upstream}--{polyA_nt}-{buffer}-{motif}--{downstream}')

# #         print(f">{location[1]}{location[2]}")
# #         print(f'{upstream}{polyA_nt}{buffer}{motif}{downstream}')



# # def extract_polyA_neighborhood_seq(location, record):

# #     # unpacking location and getting contig seq
# #     polyA_id, polyA_contig, polyA_site, polyA_strand = location
# #     contig_sequence = str(record.seq)

# #     # convert polyA_site (1-indexed) to 0-indexed
# #     polyA_site = polyA_site - 1

# #     # extract polyA neighborhood seq,
# #     # such that the polyA site is EXACTLY in the middle!
# #     ## this ensures that after reverse complementing the sequence for - strand polyA sites,
# #     ## the polyA site is still in the same position within the neighborhood as for + strand polyA sites
# #     polyA_region_seq = contig_sequence[polyA_site-24:polyA_site+25].upper()
# #     if polyA_strand == '-':
# #         polyA_region_seq = Seq(polyA_region_seq).reverse_complement()

# #     return polyA_region_seq




# # if __name__ == '__main__':
# #     main(args)





# #!/usr/bin/env python

# # standard modules
# import argparse

# # stuff you need to install
# # conda create -n <env_name> -c bioconda gffutils biopython weblogo
# import gffutils
# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord

# # argument parser
# parser = argparse.ArgumentParser(
#             description=
#         '''
#         Parses a GFF3 file for intron features (currently assumes
#         that these are written explicitly), looks up their DNA
#         sequences in the corresponding FASTA file, and outputs
#         the sequence surrounding the intron boundaries.

#         The resulting FASTA file is then used to create a
#         SequenceLogo in PNG format with the weblogo package
#         '''
#         )


# parser.add_argument(
#         "-f", "--fasta",
#         dest='fasta_file',
#         type=str,
#         required=True,
#         help="Input FASTA file")


# import weblogo as wl

# # store sequences in memory
# fin = open(test_header)
# seqs = wl.read_seq_data(fin)
# logo_data = wl.LogoData.from_seqs(seqs)

# ## setup how you want your SequenceLogo to look like
# # title = f'{intron_type}: {len(my_records)}x' 
# # logo_options = LogoOptions(show_fineprint=False, creator_text='RCL', logo_title=title)
# logo_options = wl.LogoOptions(show_fineprint=False, color_scheme=wl.classic, logo_title=title)
# ## other options include e.g. logo_title, resolution
# ## see also https://weblogo.readthedocs.io/en/latest/logo.html
# ## other colorschemes are listed in wl.std_color_schemes (monochrome, base_pairing, classic, hydrophobicity, chemistry, charge)

# # format, create and print PNG
# logo_format = wl.LogoFormat(logo_data, logo_options)
# png = wl.png_print_formatter(logo_data, logo_format)
# out_file_name = args.out_prefix + '.' + file_name.replace('fasta','png')
# with open(out_file_name, 'wb') as out:
#     out.write(png)

#!/usr/bin/env python

# Import required modules
import argparse
from Bio import SeqIO
import weblogo as wl
import io  # Import io for in-memory file handling

# Argument parser setup
parser = argparse.ArgumentParser(
)
parser.add_argument(
    "-f", "--fasta",
    dest='fasta_file',
    type=str,
    required=True,
    help="Input FASTA file containing sequences for sequence logo generation"
)
parser.add_argument(
    "-o", "--out",
    dest='out_prefix',
    type=str,
    required=True,
    help="Output prefix for the sequence logo (PNG format)"
)
args = parser.parse_args()

def main(args):
    # # Load sequences from the input FASTA file
    # fasta_records = SeqIO.parse(args.fasta_file, "fasta")
    # sequences = [str(record.seq) for record in fasta_records]

    # # Convert sequences to a single newline-separated string and use io.StringIO
    # seq_data = "\n".join(sequences)
    # seqs = wl.read_seq_data(io.StringIO(seq_data))
    fin = open(args.fasta_file)
    seqs = wl.read_seq_data(fin)

    # Create logo data
    logo_data = wl.LogoData.from_seqs(seqs)

    # Set up options for the sequence logo
    logo_options = wl.LogoOptions()
    logo_options.title = "Sequence Logo"
    logo_options.color_scheme = wl.std_color_schemes['classic']  # Example color scheme

    # Format and generate the sequence logo
    logo_format = wl.LogoFormat(logo_data, logo_options)
    png = wl.png_print_formatter(logo_data, logo_format)

    # Save the PNG file
    out_file_name = f"{args.out_prefix}.png"
    with open(out_file_name, 'wb') as out:
        out.write(png)

if __name__ == '__main__':
    main(args)

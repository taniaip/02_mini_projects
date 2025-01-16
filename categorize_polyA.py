import argparse
import gffutils
from Bio import SeqIO

def assign_polyA_to_genes(gene_pairs, polyA_sites: list[list[str, int, str, str]], db_genes, sequences):
    """
    Assign polyA sites to genes based on the criteria described, including the polyA ID in the results.
    Also adds a 'yes' or 'no' field indicating the presence of a stop codon.
    Tracks unmatched polyA sites.
    """
    results2 = []
    unmatched_polyA_ids = []  # Store unmatched polyA IDs

    for polyA in polyA_sites:
        contig, start, polyA_strand, polyA_id = polyA  # Include polyA ID here

        # Find the matching lists for this contig
        matching_gene_pairs = [gp for gp in gene_pairs if str(gp[0]) == str(contig)]
        matched = False  # Flag to track if the polyA was matched

        for gene_pair in matching_gene_pairs:
            contig, end_left, start_right, gene_left, gene_right, strand_left, strand_right = gene_pair

            # Check for the first half-interval
            if end_left == 0 and start < start_right + 3:
                if polyA_strand == "-" and strand_right == "-":
                    gene = db_genes[gene_right]
                    stop_codon = has_stop_codon(gene, sequences)
                    results2.append([contig, polyA_id, gene_right, "-", stop_codon])
                    matched = True

            elif start_right != 0:
                # For polyA and gene both on "+" strand: shift interval 3 to the left
                if polyA_strand == "+" and strand_left == "+":
                    if (end_left - 3) < start < start_right:  # Shift interval 3 to the left
                        gene = db_genes[gene_left]
                        stop_codon = has_stop_codon(gene, sequences)
                        results2.append([contig, polyA_id, gene_left, "+", stop_codon])
                        matched = True

                # For polyA and gene both on "-" strand: shift interval 3 to the right
                elif polyA_strand == "-" and strand_right == "-":
                    if end_left < start < (start_right + 3):  # Shift interval 3 to the right
                        gene = db_genes[gene_right]
                        stop_codon = has_stop_codon(gene, sequences)
                        results2.append([contig, polyA_id, gene_right, "-", stop_codon])
                        matched = True

            # Check for the last half-interval
            elif start_right == 0 and end_left - 3 < start:
                if polyA_strand == "+" and strand_left == "+":
                    gene = db_genes[gene_left]
                    stop_codon = has_stop_codon(gene, sequences)
                    results2.append([contig, polyA_id, gene_left, "+", stop_codon])
                    matched = True

        if not matched:
            unmatched_polyA_ids.append(polyA_id)  # Record unmatched polyA ID

    return results2, unmatched_polyA_ids

##ibelieve the error is rising from ere that it is not considering the reverse comliment of the negative strands
# def has_stop_codon(gene, sequences):
#     """
#     Check if a given gene has a stop codon at the end.
#     """
#     sequence = sequences[gene.seqid][gene.start - 1:gene.end].upper()
#     return "yes" if sequence[-3:] in ["TAA", "TGA", "TAG"] else "no"

def has_stop_codon(gene, sequences):
    """
    Check if a given gene has a stop codon at the end, considering the strand direction.
    """
    # Extract the gene sequence
    sequence = sequences[gene.seqid][gene.start - 1:gene.end].upper()

    if gene.strand == "+":
        # Check for stop codon on the positive strand
        return "yes" if sequence[-3:] in ["TAA", "TGA", "TAG"] else "no"
    elif gene.strand == "-":
        # Check for stop codon on the negative strand (reverse complement)
        # reverse_complement = {
        #     "TAA": "AAT",
        #     "TAG": "GAT",
        #     "TGA": "AGT"
        # }
        return "yes" if sequence[:3] in ["TTA", "TCA", "CTA"] else "no"
        
    else:
        # If strand information is unavailable, return "unknown"
        return "unknown"

def process_gff3(db):
    """
    Process the genes GFF3 database and generate gene pair information.
    """
    results = []
    genes_by_contig = {}

    for gene in db.features_of_type("gene"):
        if gene.seqid not in genes_by_contig:
            genes_by_contig[gene.seqid] = []
        genes_by_contig[gene.seqid].append(gene)

    for contig, genes in genes_by_contig.items():
        genes = sorted(genes, key=lambda g: g.start)

        first_gene = genes[0]
        results.append([contig, 0, first_gene.start, None, first_gene.id, None, first_gene.strand])

        for i in range(len(genes) - 2):
            gene_left = genes[i]
            gene_right = genes[i + 1]
            results.append([contig, gene_left.end, gene_right.start, gene_left.id, gene_right.id, gene_left.strand, gene_right.strand])

        last_gene = genes[-1]
        results.append([contig, last_gene.end, 0, last_gene.id, "", last_gene.strand, ""])

    return results


def process_polyA_sites(db):
    """
    Process the polyA GFF3 database and extract polyA site information, including the polyA ID.
    """
    polyA_sites = []
    for feature in db.features_of_type("polyA"):
        contig = feature.seqid
        start = feature.start
        strand = feature.strand
        polyA_id = feature.id
        polyA_sites.append([contig, start, strand, polyA_id])
    return polyA_sites


def filter_polyAs_by_ids(gff3_file, polyA_ids):
    """
    Extract GFF3 lines matching the given polyA IDs.
    """
    filtered_lines = []
    with open(gff3_file, "r") as gff3:
        for line in gff3:
            if line.startswith("#"):
                continue
            columns = line.strip().split("\t")
            if len(columns) > 8:
                attributes = columns[8]
                polyA_id = None
                for attr in attributes.split(";"):
                    if attr.startswith("ID="):
                        polyA_id = attr.split("=")[1]
                        break
                if polyA_id and polyA_id in polyA_ids:
                    filtered_lines.append(line)
    return filtered_lines


def load_fasta_sequences(fasta_file):
    """
    Load sequences from a FASTA file into a dictionary.
    """
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


def main():
    parser = argparse.ArgumentParser(description="Categorize polyA sites into 'yes', 'no', and 'not matched'.")
    parser.add_argument("-g1", "--genes", required=True, help="Input GFF3 file with gene information.")
    parser.add_argument("-g2", "--polyA", required=True, help="Input GFF3 file with polyA site information.")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file containing genomic sequences.")
    parser.add_argument("--output_yes", required=True, help="Output GFF3 file for 'yes' polyA entries.")
    parser.add_argument("--output_no", required=True, help="Output GFF3 file for 'no' polyA entries.")
    parser.add_argument("--output_not_matched", required=True, help="Output GFF3 file for unmatched polyA entries.")
    args = parser.parse_args()

    db_genes = gffutils.create_db(args.genes, dbfn=":memory:", merge_strategy="create_unique", keep_order=True)
    db_polyA = gffutils.create_db(args.polyA, dbfn=":memory:", merge_strategy="create_unique", keep_order=True)
    sequences = load_fasta_sequences(args.fasta)

    gene_pairs = process_gff3(db_genes)
    polyA_sites = process_polyA_sites(db_polyA)

    categorized_polyA, unmatched_polyA_ids = assign_polyA_to_genes(gene_pairs, polyA_sites, db_genes, sequences)

    yes_ids = {entry[1] for entry in categorized_polyA if entry[-1] == "yes"}
    no_ids = {entry[1] for entry in categorized_polyA if entry[-1] == "no"}

    yes_lines = filter_polyAs_by_ids(args.polyA, yes_ids)
    no_lines = filter_polyAs_by_ids(args.polyA, no_ids)
    unmatched_lines = filter_polyAs_by_ids(args.polyA, unmatched_polyA_ids)

    with open(args.output_yes, "w") as yes_file, open(args.output_no, "w") as no_file, open(args.output_not_matched, "w") as unmatched_file:
        yes_file.writelines(yes_lines)
        no_file.writelines(no_lines)
        unmatched_file.writelines(unmatched_lines)


if __name__ == "__main__":
    main()




##ATTENTION!!!! It is working perfectly for yes and no groups(not considering non-matched!):

# import argparse
# import gffutils
# import pprint
# from Bio import SeqIO

# # Custom PrettyPrinter with increased width
# pprint.pprint = lambda obj, **kwargs: pprint.PrettyPrinter(width=120).pprint(obj)

# parser = argparse.ArgumentParser()
# parser.add_argument(
#     "-g1", "--genes",
#     dest="genes_gff3_file",
#     type=str,
#     required=True,
#     help="Input GFF3 file with gene information."
# )
# parser.add_argument(
#     "-g2", "--polyA",
#     dest="polyA_gff3_file",
#     type=str,
#     required=True,
#     help="Input GFF3 file with polyA site information."
# )
# parser.add_argument(
#     "-f", "--fasta",
#     dest="fasta_file",
#     type=str,
#     required=True,
#     help="Input FASTA file containing genomic sequences."
# )
# parser.add_argument(
#     "--output_yes",
#     type=str,
#     required=True,
#     help="Output GFF3 file for 'yes' polyA entries."
# )
# parser.add_argument(
#     "--output_no",
#     type=str,
#     required=True,
#     help="Output GFF3 file for 'no' polyA entries."
# )


# def load_fasta_sequences(fasta_file):
#     """
#     Load sequences from a FASTA file into a dictionary.
#     """
#     sequences = {}
#     for record in SeqIO.parse(fasta_file, "fasta"):
#         sequences[record.id] = str(record.seq)
#     return sequences


# def process_gff3(db):
#     """
#     Process the genes GFF3 database and generate gene pair information.
#     """
#     results = []
#     genes_by_contig = {}

#     # Group genes by contig
#     for gene in db.features_of_type("gene"):
#         if gene.seqid not in genes_by_contig:
#             genes_by_contig[gene.seqid] = []  # Initialize list for this contig
#         genes_by_contig[gene.seqid].append(gene)

#     # Process each contig
#     for contig, genes in genes_by_contig.items():
#         genes = sorted(genes, key=lambda g: g.start)  # Sort genes by start position

#         # First gene
#         first_gene = genes[0]
#         results.append([
#             contig,
#             0,
#             first_gene.start,
#             None,
#             first_gene.id,
#             None,
#             first_gene.strand,
#         ])

#         # All genes in between first and last gene
#         for i in range(len(genes) - 2):
#             gene_left = genes[i]
#             gene_right = genes[i + 1]

#             results.append([
#                 contig,
#                 gene_left.end,
#                 gene_right.start,
#                 gene_left.id,
#                 gene_right.id,
#                 gene_left.strand,
#                 gene_right.strand
#             ])

#         # Last gene
#         last_gene = genes[-1]
#         results.append([
#             contig,
#             last_gene.end,
#             0,
#             last_gene.id,
#             "",
#             last_gene.strand,
#             ""
#         ])

#     return results


# def process_polyA_sites(db) -> list[list[str, int, str, str]]:
#     """
#     Process the polyA GFF3 database and extract polyA site information, including the polyA ID.
#     """
#     polyA_sites = []
#     for feature in db.features_of_type("polyA"):
#         contig = feature.seqid
#         start = feature.start
#         strand = feature.strand
#         polyA_id = feature.id  # Extract the polyA ID
#         polyA_sites.append([contig, start, strand, polyA_id])
#     return polyA_sites


# def has_stop_codon(gene, sequences):
#     """
#     Check if a given gene has a stop codon at the end.
#     """
#     # Extract the gene sequence
#     sequence = sequences[gene.seqid][gene.start - 1:gene.end] 

#     # Convert to uppercase for consistency
#     sequence = sequence.upper()

#     # Check the last three bases for stop codons
#     return "yes" if sequence[-3:] in ["TAA", "TGA", "TAG"] else "no"


# def assign_polyA_to_genes(gene_pairs, polyA_sites: list[list[str, int, str, str]], db_genes, sequences):
#     """
#     Assign polyA sites to genes based on the criteria described, including the polyA ID in the results.
#     Also adds a 'yes' or 'no' field indicating the presence of a stop codon.
#     """
#     results2 = []

#     for polyA in polyA_sites:
#         contig, start, polyA_strand, polyA_id = polyA  # Include polyA ID here

#         # Find the matching lists for this contig
#         matching_gene_pairs = [gp for gp in gene_pairs if str(gp[0]) == str(contig)]

#         for gene_pair in matching_gene_pairs:
#             contig, end_left, start_right, gene_left, gene_right, strand_left, strand_right = gene_pair

#             # Check for the first half-interval
#             if end_left == 0 and start < start_right + 3:
#                 if polyA_strand == "-" and strand_right == "-":
#                     gene = db_genes[gene_right]
#                     stop_codon = has_stop_codon(gene, sequences)
#                     results2.append([contig, polyA_id, gene_right, "-", stop_codon])

#             elif start_right != 0:
#                 # For polyA and gene both on "+" strand: shift interval 3 to the left
#                 if polyA_strand == "+" and strand_left == "+":
#                     if (end_left - 3) < start < start_right:  # Shift interval 3 to the left
#                         gene = db_genes[gene_left]
#                         stop_codon = has_stop_codon(gene, sequences)
#                         results2.append([contig, polyA_id, gene_left, "+", stop_codon])

#                 # For polyA and gene both on "-" strand: shift interval 3 to the right
#                 elif polyA_strand == "-" and strand_right == "-":
#                     if end_left < start < (start_right + 3):  # Shift interval 3 to the right
#                         gene = db_genes[gene_right]
#                         stop_codon = has_stop_codon(gene, sequences)
#                         results2.append([contig, polyA_id, gene_right, "-", stop_codon])

#             # Check for the last half-interval
#             elif start_right == 0 and end_left - 3 < start:
#                 if polyA_strand == "+" and strand_left == "+":
#                     gene = db_genes[gene_left]
#                     stop_codon = has_stop_codon(gene, sequences)
#                     results2.append([contig, polyA_id, gene_left, "+", stop_codon])

#     return results2


# def filter_polyAs_by_ids(gff3_file, gene_polyA_list, output_yes, output_no):
#     """
#     Write polyA entries from the GFF3 file into separate 'yes' and 'no' output GFF3 files
#     based on their IDs in the gene_polyA_list.
#     """
#     yes_ids = {entry[1] for entry in gene_polyA_list if entry[-1] == "yes"}
#     no_ids = {entry[1] for entry in gene_polyA_list if entry[-1] == "no"}

#     with open(gff3_file, "r") as gff3, open(output_yes, "w") as yes_file, open(output_no, "w") as no_file:
#         for line in gff3:
#             if line.startswith("#"):  # Write header lines
#                 yes_file.write(line)
#                 no_file.write(line)
#                 continue

#             columns = line.strip().split("\t")
#             if len(columns) > 8:  # Ensure the line has attributes
#                 attributes = columns[8]
#                 polyA_id = None
#                 for attr in attributes.split(";"):
#                     if attr.startswith("ID="):  # Extract the ID attribute
#                         polyA_id = attr.split("=")[1]
#                         break
#                 if polyA_id in yes_ids:
#                     yes_file.write(line)
#                 elif polyA_id in no_ids:
#                     no_file.write(line)


# def main(args):
#     # Create a database for the genes GFF3 file
#     db_genes = gffutils.create_db(
#         args.genes_gff3_file,
#         dbfn=":memory:",
#         merge_strategy="create_unique",
#         keep_order=True,
#     )
#     # Create a database for the polyA sites GFF3 file
#     db_polyA = gffutils.create_db(
#         args.polyA_gff3_file,
#         dbfn=":memory:",
#         merge_strategy="create_unique",
#         keep_order=True,
#     )

#     # Load sequences from FASTA
#     sequences = load_fasta_sequences(args.fasta_file)

#     # Generate gene pairs from the genes GFF3 file
#     result_list_gene_pairs = process_gff3(db_genes)

#     # Extract polyA sites from the polyA sites GFF3 file
#     polyA_sites = process_polyA_sites(db_polyA)

#     # Assign polyA sites to genes
#     gene_polyA_list = assign_polyA_to_genes(result_list_gene_pairs, polyA_sites, db_genes, sequences)

#     # Filter and write 'yes' and 'no' polyA entries to separate GFF3 files
#     filter_polyAs_by_ids(args.polyA_gff3_file, gene_polyA_list, args.output_yes, args.output_no)


# if __name__ == "__main__":
#     args = parser.parse_args()
#     main(args)








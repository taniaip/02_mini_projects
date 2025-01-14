import argparse
import gffutils

import pprint

# Custom PrettyPrinter with increased width(chatGPT)
pprint.pprint = lambda obj, **kwargs: pprint.PrettyPrinter(width=120).pprint(obj)


parser = argparse.ArgumentParser()
parser.add_argument(
    "-g", "--gff3",
    dest="gff3_file",
    type=str,
    required=True,
    help="Input GFF3 file with gene information."
)


def process_gff3(db):
    results = []
    genes_by_contig = {}

    # Group genes by contig
    for gene in db.features_of_type("gene"):
        if gene.seqid not in genes_by_contig:
            genes_by_contig[gene.seqid] = []  # Initialize list for this contig
        genes_by_contig[gene.seqid].append(gene)

    # Process each contig
    for contig, genes in genes_by_contig.items():
        genes = sorted(genes, key=lambda g: g.start)  # Sort genes by start position

        first_gene = genes[0]
        results.append([
            contig,
            0,
            first_gene.start,
            None,  
            first_gene.id,
            None,  
            first_gene.strand,
        ])

        for i in range(len(genes) - 2):
            gene_left = genes[i]
            gene_right = genes[i + 1]

            results.append([
                contig,
                gene_left.end,
                gene_right.start,
                gene_left.id,
                gene_right.id,
                gene_left.strand,
                gene_right.strand
            ])
        last_gene = genes[-1]
        results.append([
            contig,
            last_gene.end,
            0,
            last_gene.id,
            "",
            last_gene.strand,
            ""
            ])

    return results 


def main(args):
    # Create a database from the GFF3 file
    db = gffutils.create_db(
        args.gff3_file,
        dbfn=":memory:",
        merge_strategy="create_unique"
    )

    # Process the database and generate results
    result_list = process_gff3(db)

    # Pretty-print the results
    pprint.pprint(result_list)


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)






# """
# combind code:
# """

# import argparse
# from pprint import pprint
# import gffutils

# parser = argparse.ArgumentParser(description="Match polyA sites to genes based on GFF3 data.")
# parser.add_argument(
#     "-g", "--gff3",
#     dest="gff3_file",
#     type=str,
#     required=True,
#     help="Input GFF3 file with gene information and polyA sites."
# )

# def process_gene_pairs(db):
#     """
#     Process the GFF3 database and extract gene pairs by contig.
#     """
#     results = []
#     genes_by_contig = {}

#     # Group genes by contig
#     for gene in db.features_of_type("gene"):
#         genes_by_contig.setdefault(gene.seqid, []).append(gene)

#     # Process each contig
#     for contig, genes in genes_by_contig.items():
#         # Sort genes by start position
#         genes = sorted(genes, key=lambda g: g.start)
        
#         first_gene = genes[0]
#         results.append([
#             contig,
#             0,
#             first_gene_gene.start,
#             None,  
#             first_gene.id,
#             None,  
#             first_gene.strand,
#         ])
#         # Pair genes
#         for i in range(len(genes) - 1):
#             gene_left = genes[i]
#             gene_right = genes[i + 1]
#             results.append([
#                 contig,
#                 gene_left.end,
#                 gene_right.start,
#                 gene_left.id,
#                 gene_right.id,
#                 gene_left.strand,
#                 gene_right.strand,
#             ])


#         last_gene = genes[-1]
#         results.append([
#             contig,
#             last_gene.end,
#             0,  # No next gene start
#             last_gene.id,
#             None,  # No right gene ID
#             last_gene.strand,
#             None,  # No right gene strand
#         ])
    
#     return results

# def process_polyA_sites(db):
#     """
#     Process the GFF3 database and extract polyA site information.
#     """
#     return [
#         [feature.seqid, feature.end, feature.strand]
#         for feature in db.features_of_type("polyA")
#     ]

# def match_polyA_to_genes(gene_pairs, polyA_sites):
#     """
#     Match polyA sites to genes based on contig, position, and strand.
#     """
#     results = []

#     for polyA in polyA_sites:
#         contig, start, polyA_strand = polyA

#         # Find gene pairs for the current contig
#         matching_gene_pairs = [gp for gp in gene_pairs if str(gp[0]) == str(contig)]

#         for gene_pair in matching_gene_pairs:
#             _, end_left, start_right, gene_left, gene_right, strand_left, strand_right = gene_pair

#             # Match polyA site to a gene
#             if start_right != 0 and end_left < start < start_right:
#                 if polyA_strand == "+" and strand_left == "+":
#                     results.append([contig, gene_left, start, "+"])
#                 elif polyA_strand == "-" and strand_right == "-":
#                     results.append([contig, gene_right, start, "-"])
#                 break
#             elif start_right == 0 and end_left < start:  # Last gene case
#                 if polyA_strand == "+" and strand_left == "+":
#                     results.append([contig, gene_left, start, "+"])
    
#     return results

# def main(args):
#     # Create a GFF3 database
#     db = gffutils.create_db(
#         args.gff3_file,
#         dbfn=":memory:",
#         merge_strategy="create_unique"
#     )

#     # Process gene pairs and polyA sites
#     gene_pairs = process_gene_pairs(db)
#     polyA_sites = process_polyA_sites(db)

#     # Match polyA sites to genes
#     results = match_polyA_to_genes(gene_pairs, polyA_sites)

#     # Output results
#     pprint(results)

# if __name__ == "__main__":
#     args = parser.parse_args()
#     main(args)



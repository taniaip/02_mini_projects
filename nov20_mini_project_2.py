import argparse
from pprint import pprint
import gffutils
from Bio import SeqIO
from Bio.Seq import Seq


parser = argparse.ArgumentParser()
parser.add_argument(
    "-g", "--gff3",
    dest="gff3_file",
    type=str,
    required=True,
    help="Input GFF3 file with polyA sites"
)
parser.add_argument(
    "-gp", "--gene_pairs_file",
    dest="gene_pairs_file",
    type=str,
    required=True,
    help="Input file containing gene pairs"
)


def read_gene_pairs(file_path):
    """
    Read gene pairs from the text file into a list.
    """
    gene_pairs = []
    with open(file_path, "r") as f:
        gene_pairs = [
            [parts[0], int(parts[1]), int(parts[2]),
             parts[3], parts[4], parts[5], parts[6]]
            for line in f if (parts := line.strip().split("\t"))
        ]
    return gene_pairs

def format_gene_pairs_with_newlines(gene_pairs):
    """
    Format gene pairs with newline-separated rows.
    """
    formatted_pairs = "\n".join("\t".join(map(str, row)) for row in gene_pairs)
    return formatted_pairs + "\n"


def assign_polyA_to_genes(gene_pairs, polyA_sites: list[list[str,int,str]]):
    """
    Assign polyA sites to genes based on the criteria described.
    """
    results = []

    for polyA in polyA_sites:
        contig, start, polyA_strand = polyA

        # Find the matching lists for this contig
        matching_gene_pairs = [gp for gp in gene_pairs if str(gp[0]) == str(contig)]

        for gene_pair in matching_gene_pairs:
            _, end_left, start_right, gene_left, gene_right, strand_left, strand_right = gene_pair
            #check for the first half-interval    
            if end_left == 0 and start < start_right:
                if polyA_strand == "-" and strand_right == "-":
                    results.append([contig, gene_right, "-"])
                else:
                    break                 

            # Check if polyA start site falls into the interval
            elif start_right != 0 and end_left < start < start_right:
                if polyA_strand == "+" and strand_left == "+":
                    results.append([contig, gene_left, "+"])
                elif polyA_strand == "-" and strand_right == "-":
                    results.append([contig, gene_right, "-"])
                else:
                    break   

            #check for the last half-interval         
            elif start_right == 0 and end_left < start:
                if polyA_strand == "+" and strand_left == "+":
                    results.append([contig, gene_left, "+"])
                else:
                    break  
            else:
                break          

    return results


def process_polyA_sites(db) -> list[list[str,int,str]]:
    """
    Process the GFF3 database and extract polyA site information.
    """
    polyA_sites = []
    for feature in db.features_of_type("polyA"):
        contig = feature.seqid
        start = feature.start
        strand = feature.strand
        polyA_sites.append([contig, start, strand])
    return polyA_sites


def main(args):
    db = gffutils.create_db(
        args.gff3_file,
        dbfn=":memory:",
        merge_strategy="create_unique"
    )

    polyA_sites: list[list[str,int,str]] = process_polyA_sites(db)

    gene_pairs = read_gene_pairs(args.gene_pairs_file)

    results = assign_polyA_to_genes(gene_pairs, polyA_sites)

    pprint(results)


if __name__ == "__main__":
    args = parser.parse_args()
    main(args)

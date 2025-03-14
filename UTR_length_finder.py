# $ python UTR_length_finder.py -g1 ST2_orfs.gff3 -g2 march4_has_stop_polyA.gff3 -f ST2_sorted_masked.fasta > UTR_length_info_new.txt

# $ grep -c "there must be" UTR_length_info_new.txt
# 52 


# $ python UTR_length_finder.py -g1 ST2_orfs.gff3 -g2 march4_has_stop_polyA.gff3 -f ST2_sorted_masked.fasta > UTR_length_info_new.txt
#!/usr/bin/env python
import sys
import argparse

import gffutils
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


from typing import List, Dict, Optional
from gffutils.feature import Feature
from gffutils.interface import FeatureDB
from Bio.SeqRecord import SeqRecord


# Argument Parser
parser = argparse.ArgumentParser(description="Categorize polyA sites into genomic contexts.")
parser.add_argument("-g1", "--genes", required=True, help="Input GFF3 file with gene information.")
parser.add_argument("-g2", "--polyA", required=True, help="Input GFF3 file with polyA sites which are associated with genes having stop codons.")
parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file containing genomic sequences.")
args = parser.parse_args()


class Interval:
    """Class representing an interval between two genes."""
    def __init__(self, contig: str, gene_left: Optional[Feature], gene_right: Optional[Feature]):
        self.contig = contig
        self.start = gene_left.end if gene_left else 0
        self.end = gene_right.start if gene_right else sys.maxsize
        self.gene_left = gene_left
        self.gene_right = gene_right
        self.strand_left = gene_left.strand if gene_left else None
        self.strand_right = gene_right.strand if gene_right else None


def main() -> None:
    # FIX: Use in-memory DB cautiously, large files might cause issues.
    db_genes: FeatureDB = gffutils.create_db(args.genes, dbfn="genes.db", force=True, merge_strategy="create_unique")
    db_polyA: FeatureDB = gffutils.create_db(args.polyA, dbfn="polyA.db", force=True, merge_strategy="create_unique")
    sequences = SeqIO.index(args.fasta, 'fasta')

    # Generate intervals
    gene_intervals = process_gene_intervals(db_genes)

    UTR_lengths_list = find_UTR_lengths(gene_intervals, db_polyA, sequences)

    df = pd.DataFrame(UTR_lengths_list, columns=["UTR Length"])
    plt.xlim(-10, 200)  # Set max x-axis to 1000
    sns.histplot(df, x="UTR Length", bins=1000, kde=True).get_figure().savefig("UTR_length_distribution.jpg", dpi=300)
      

def process_gene_intervals(db: FeatureDB) -> List[Interval]:
    """
    Generate intervals between 'big' genes in the GFF3 database,
    skipping smaller genes overlapped by a bigger gene.
    """
    intervals: List[Interval] = []
    genes_by_contig: Dict[str, List[Feature]] = {}

    # Collect genes by contig
    for gene in db.features_of_type("gene"):
        if gene.seqid not in genes_by_contig:
            genes_by_contig[gene.seqid] = []
        genes_by_contig[gene.seqid].append(gene)

    # Process each contig individually
    for contig, genes in genes_by_contig.items():
        # Sort genes by their start positions and in case of starting at the same point prioritise the ones with further ends
        genes.sort(key=lambda g: (g.start, -g.end))

        # Create interval from "start of genome" to the first gene
        first_gene = genes[0]
        intervals.append(Interval(contig=contig, gene_left=None, gene_right=first_gene))
        correct_start = first_gene

        # Iterate over the remaining genes
        for i in range(1, len(genes)):
            # If next gene is overshadowed by correct_start (start < correct_start.end), skip it
            if genes[i].start < correct_start.end:
                continue

            # Otherwise, create an interval from the last big gene to this new gene
            intervals.append(Interval(contig=contig, gene_left=correct_start, gene_right=genes[i]))

            correct_start = genes[i]

        # Create interval from the last big gene to "end of genome"
        intervals.append(Interval(contig=contig, gene_left=correct_start, gene_right=None))

    return intervals


def has_stop_codon(gene: Feature, sequences: Dict[str, SeqRecord]) -> bool:
    """Check if a gene has a stop codon based on its strand."""
    sequence = str(sequences[gene.seqid].seq[gene.start - 1: gene.end]).upper()
    if gene.strand == "+":
        return sequence[-3:] in ["TAA", "TGA", "TAG"]
    if gene.strand == "-":
        return sequence[:3] in ["TTA", "TCA", "CTA"]
    return False


def find_UTR_lengths(
    gene_intervals: List[Interval],
    db_polyA: FeatureDB,
    sequences: Dict[str, SeqRecord],
) -> List[int]: 
    """Assign polyA sites to genes, categorize them, and return UTR lengths."""
    UTR_lengths = []
    geneIdToCheck = None
    UTRlengthToCheck = None
    for polyA in db_polyA.features_of_type('polyA'):
        for interval in gene_intervals:
            if polyA.seqid != interval.contig:
                continue

            # Check if polyA falls within an interval and the strands match
            if interval.start - 3 <= polyA.start <= interval.end + 3:
                if polyA.strand == "-" and interval.strand_right == "-":
                    current_gene = interval.gene_right.id
                    current_UTR_length = interval.end - polyA.start


                elif polyA.strand == "+" and interval.strand_left == "+":
                    current_gene = interval.gene_left.id
                    current_UTR_length =polyA.start - interval.start


                if current_UTR_length < 1000:
                    UTR_lengths.append(current_UTR_length)
                else:
                    print(f"WARNING: UTR ({current_UTR_length} bp) for gene {current_gene} is unusually long; a missing gene may be possible.")

                print(current_gene, polyA.id, current_UTR_length)

                # Only do the comparison if both geneIdToCheck and UTRlengthToCheck are already defined
                if geneIdToCheck is not None and UTRlengthToCheck is not None:
                # Check if we have the same gene, and if the difference is bigger than 50
                    if geneIdToCheck == current_gene and abs(UTRlengthToCheck - current_UTR_length) > 50:
                        print(f"WARNING: there might be a gene missing between {polyAtoCheck} and {polyA.id}")

                # Update these variables for the next comparison
                geneIdToCheck = current_gene
                UTRlengthToCheck = current_UTR_length
                polyAtoCheck = polyA.id

    return UTR_lengths


if __name__ == "__main__":
    main()

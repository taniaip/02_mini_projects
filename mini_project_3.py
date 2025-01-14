import argparse
import gffutils
from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument(
        "-g", "--gff3",
        dest='gff3_file',
        type=str,
        required=True,
        help="Input GFF3 genome file")

# Parse the arguments
args = parser.parse_args()

# load a GFF file and create an sqlite3 database file
# db = gffutils.create_db(args.gff3_file, 'file.db')

# load a GFF file and store the sqlite3 database into memory
db = gffutils.create_db(args.gff3_file, ':memory:')
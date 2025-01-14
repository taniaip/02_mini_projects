import argparse


def load_polyA_ids(file_path):
    """
    Load polyA IDs from a file (e.g., 'yes' or 'no' text file).
    Extract the second element (polyA ID) from each line.
    """
    polyA_ids = set()
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("[[") or not line.startswith("["):
                continue  # Skip invalid or improperly formatted lines
            line = line.lstrip("[[").rstrip("]],").strip()  # Strip surrounding brackets
            elements = line.split(",")  # Split by commas
            if len(elements) > 1:
                polyA_id = elements[1].strip().strip("'")  # Extract the second element
                polyA_ids.add(polyA_id)
    return polyA_ids


def filter_polyAs_by_ids(gff3_file, polyA_ids, output_file):
    """
    Write polyA entries from the GFF3 file to the output file
    if their ID matches the given polyA IDs.
    """
    with open(gff3_file, "r") as gff3, open(output_file, "w") as output:
        for line in gff3:
            if line.startswith("#"):  # Write header lines
                output.write(line)
                continue
            columns = line.strip().split("\t")
            if len(columns) > 8:  # Ensure the line has attributes
                attributes = columns[8]
                polyA_id = None
                for attr in attributes.split(";"):
                    if attr.startswith("ID="):  # Extract the ID attribute
                        polyA_id = attr.split("=")[1]
                        break
                if polyA_id and polyA_id in polyA_ids:
                    output.write(line)  # Write the matched polyA line to the output file


def main():
    parser = argparse.ArgumentParser(description="Generate separate GFF3 files for 'yes' and 'no' polyA IDs.")
    parser.add_argument("-g", "--gff3", required=True, help="Input GFF3 file with polyA site information.")
    parser.add_argument("-y", "--yes", required=True, help="File containing 'yes' polyA IDs.")
    parser.add_argument("-n", "--no", required=True, help="File containing 'no' polyA IDs.")
    parser.add_argument("--output_yes", required=True, help="Output GFF3 file for 'yes' polyA entries.")
    parser.add_argument("--output_no", required=True, help="Output GFF3 file for 'no' polyA entries.")
    args = parser.parse_args()

    # Load polyA IDs from 'yes' and 'no' files
    yes_ids = load_polyA_ids(args.yes)
    no_ids = load_polyA_ids(args.no)

    # Generate GFF3 files for 'yes' and 'no'
    filter_polyAs_by_ids(args.gff3, yes_ids, args.output_yes)
    filter_polyAs_by_ids(args.gff3, no_ids, args.output_no)


if __name__ == "__main__":
    main()

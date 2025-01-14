import argparse


def load_polyA_ids_from_text(file_path):
    """
    Extract polyA IDs (second element in the list) from the text file.
    """
    polyA_ids = set()
    with open(file_path, "r") as f:
        for line in f:
            # Process each line manually without literal_eval
            line = line.strip()
            if line.startswith("[[") or not line.startswith("["):
                continue  # Skip invalid or improperly formatted lines
            line = line.lstrip("[[").rstrip("]],").strip()  # Strip surrounding brackets
            elements = line.split(",")  # Split by commas
            if len(elements) > 1:
                polyA_id = elements[1].strip().strip("'")  # Extract the second element
                polyA_ids.add(polyA_id)
    return polyA_ids


def filter_unmatched_polyAs(gff3_file, matched_polyA_ids, output_file):
    """
    Write unmatched polyA entries from the GFF3 file to the output file.
    """
    with open(gff3_file, "r") as gff3, open(output_file, "w") as output:
        for line in gff3:
            if line.startswith("#"):  # Keep header lines
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
                if polyA_id and polyA_id not in matched_polyA_ids:
                    output.write(line)  # Write the unmatched polyA line to the output file


def main():
    parser = argparse.ArgumentParser(description="Find unmatched polyA entries in GFF3 file.")
    parser.add_argument(
        "-g", "--gff3", required=True, help="Input GFF3 file with polyA site information."
    )
    parser.add_argument(
        "-t", "--text", required=True, help="Input text file with matched polyA IDs."
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output file to store unmatched polyA entries."
    )
    args = parser.parse_args()

    # Load matched polyA IDs from the text file
    matched_polyA_ids = load_polyA_ids_from_text(args.text)

    # Filter unmatched polyAs and write to output
    filter_unmatched_polyAs(args.gff3, matched_polyA_ids, args.output)


if __name__ == "__main__":
    main()

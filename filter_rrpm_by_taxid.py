import csv
import sys
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rrpm", required=True)
    parser.add_argument("--taxids", required=True)
    args = parser.parse_args()

    # Open the taxids file
    with open(args.taxids) as taxids_file:
        reader = csv.DictReader(taxids_file)

        # Create a set containing the taxids
        taxids = set(row["taxID"] for row in reader)

    # List for storing the output filtered data
    output = []

    # Open the rrpm and filter using the taxids
    with open(args.rrpm) as rrpm_file:
        reader = csv.DictReader(rrpm_file)

        # For each row in the rrpm, check if the taxid is in the set of taxids
        for row in reader:
            # If the taxid is in the set of taxids, add it to the output data
            if row["taxID"] in taxids:
                output.append(row)

    # Write the data to stdout (where it can be directed into a file)
    if output:
        writer = csv.DictWriter(sys.stdout, fieldnames=output[0])
        writer.writeheader()
        writer.writerows(output)


if __name__ == "__main__":
    main()

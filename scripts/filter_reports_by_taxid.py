import os
import csv
import sys
import glob
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--folder", required=True)
    parser.add_argument("--taxids", required=True)
    args = parser.parse_args()

    # Open the taxids file
    with open(args.taxids) as taxids_file:
        reader = csv.DictReader(taxids_file)

        # Create a set containing the taxids
        taxids = set(row["taxID"] for row in reader)

    # List of reports in the folder
    report_paths = glob.glob(os.path.join(args.folder, "*_species-level-report.tsv"))

    # List for storing the output filtered data
    output = []

    # Open each report and filter each report using the taxids
    for report_path in report_paths:
        sample_name = os.path.basename(report_path).split("_")[0]

        # Open the current report
        with open(report_path) as report_file:
            reader = csv.DictReader(report_file, delimiter="\t")

            # For each row in the report, check if the taxid is in the set of taxids
            for row in reader:
                # If the taxid is in the set of taxids, add it to the output data
                if row["taxID"] in taxids:
                    out_row = {"sampleName": sample_name}
                    out_row.update(row)
                    output.append(out_row)

    # Write the data to stdout (where it can be directed into a file)
    if output:
        writer = csv.DictWriter(sys.stdout, fieldnames=output[0])
        writer.writeheader()
        writer.writerows(output)


if __name__ == "__main__":
    main()

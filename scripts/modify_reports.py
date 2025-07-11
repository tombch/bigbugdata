import csv
import os
from pathlib import Path
import argparse


def create_complete_reports(
    species_reports: list[str],
    dna_totalreads: str,
    rna_totalreads: str,
    output_dir: str,
):
    """
    Create complete Kraken reports with simulated root and unclassified entries.
    """

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Load total reads data
    total_reads_data = {}

    # Load DNA total reads
    with open(dna_totalreads) as dna_file:
        dna_data = csv.reader(dna_file, delimiter="\t")
        for row in dna_data:
            sample = row[0]
            total_reads = int(row[2])
            sample_name = sample.rpartition("_")[0]
            total_reads_data[sample_name] = total_reads

    # Load RNA total reads (this will update/override DNA values if same sample names)
    with open(rna_totalreads) as rna_file:
        rna_data = csv.reader(rna_file, delimiter="\t")
        for row in rna_data:
            sample = row[0]
            total_reads = int(row[2])
            sample_name = sample.rpartition("_")[0]
            total_reads_data[sample_name] = total_reads

    # Process each species-level report
    for species_report in species_reports:
        # Get sample name from filename
        sample_name = os.path.basename(species_report).rpartition("_")[0]

        if sample_name not in total_reads_data:
            print(f"Warning: No total reads data found for sample {sample_name}")
            continue

        total_reads = total_reads_data[sample_name]

        # Read the species-level report and calculate classified reads
        classified_reads = 0
        species_entries = []

        with open(species_report) as report_file:
            reader = csv.DictReader(report_file, delimiter="\t")
            fieldnames = (
                reader.fieldnames
                if reader.fieldnames
                else [
                    "%",
                    "reads",
                    "taxReads",
                    "kmers",
                    "dup",
                    "cov",
                    "taxID",
                    "rank",
                    "taxName",
                ]
            )

            for row in reader:
                if row["rank"] == "species":
                    classified_reads += int(row["reads"])
                    species_entries.append(row)

        # Calculate unclassified reads
        unclassified_reads = total_reads - classified_reads

        if unclassified_reads < 0:
            print(
                f"Warning: Sample {sample_name} has more classified reads than total reads!"
            )
            unclassified_reads = 0

        # Create complete report filename
        complete_report_path = (
            Path(output_dir) / f"{sample_name}_species-level-report.tsv"
        )

        # Write complete report
        with open(complete_report_path, "w") as complete_file:
            writer = csv.DictWriter(
                complete_file, fieldnames=fieldnames, delimiter="\t"
            )
            writer.writeheader()

            # Write unclassified entry (taxID 0)
            if unclassified_reads > 0:
                unclassified_entry = {
                    "%": f"{(unclassified_reads / total_reads) * 100:.4f}",
                    "reads": str(unclassified_reads),
                    "taxReads": str(unclassified_reads),
                    "kmers": "0",
                    "dup": "0",
                    "cov": "0",
                    "taxID": "0",
                    "rank": "unclassified",
                    "taxName": "unclassified",
                }
                writer.writerow(unclassified_entry)

            # Write root entry (taxID 1)
            root_entry = {
                "%": f"{(classified_reads / total_reads) * 100:.4f}",
                "reads": str(classified_reads),
                "taxReads": str(classified_reads),
                "kmers": str(sum(int(entry["kmers"]) for entry in species_entries)),
                "dup": "0",  # Not meaningful at root level
                "cov": "0",  # Not meaningful at root level
                "taxID": "1",
                "rank": "root",
                "taxName": "root",
            }
            writer.writerow(root_entry)

            # Write all species entries
            for entry in species_entries:
                writer.writerow(entry)

        print(f"Created complete report: {complete_report_path}")
        print(f"  Total reads: {total_reads:,}")
        print(f"  Classified reads: {classified_reads:,}")
        print(f"  Unclassified reads: {unclassified_reads:,}")
        print()


def main():
    parser = argparse.ArgumentParser(
        description="Generate complete Kraken reports with simulated root and unclassified entries"
    )
    parser.add_argument(
        "--species-reports",
        required=True,
        nargs="+",
        help="Species-level Kraken report files",
    )
    parser.add_argument(
        "--dna-totalreads", required=True, help="DNA total reads TSV file"
    )
    parser.add_argument(
        "--rna-totalreads", required=True, help="RNA total reads TSV file"
    )
    parser.add_argument(
        "--output-dir",
        default="complete_reports",
        help="Output directory for complete reports (default: complete_reports)",
    )

    args = parser.parse_args()

    create_complete_reports(
        species_reports=args.species_reports,
        dna_totalreads=args.dna_totalreads,
        rna_totalreads=args.rna_totalreads,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()

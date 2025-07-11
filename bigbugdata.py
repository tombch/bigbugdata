from scipy import stats
from pathlib import Path
from typing import Any
import argparse
import csv
import os
import re


def get_output_paths(results_dir: str) -> tuple[Path, Path, Path]:
    """
    Get the output file paths for combined species, rrpm, and tophits.
    """

    # Create the results directory if it doesn't exist
    Path(results_dir).mkdir(parents=True, exist_ok=True)

    # Output file paths
    combined_species_out = Path(results_dir) / "combined_species.csv"
    rrpm_out = Path(results_dir) / "rrpm.csv"
    tophits_out = Path(results_dir) / "tophits.csv"

    return combined_species_out, rrpm_out, tophits_out


def get_sample_names(report_paths: list[str]) -> dict[str, str]:
    """
    Get mapping of sample names to their report paths by partitioning on the last underscore.
    """

    return {
        os.path.basename(report).rpartition("_")[0]: report for report in report_paths
    }


def get_negative_groups(
    sample_names: list[str],
    negative_group_patterns: list[tuple[str, str]] | None,
) -> dict[str, set[str]]:
    """
    Get a dictionary mapping negative samples to their negative groups based on the provided patterns.
    """

    negative_groups = {}

    if negative_group_patterns is not None:
        for negative_sample_pattern, negative_group_pattern in negative_group_patterns:
            # Find sample that matches the negative sample pattern
            matching_negative_samples = [
                sample
                for sample in sample_names
                if re.search(negative_sample_pattern, sample)
            ]
            if len(matching_negative_samples) != 1:
                raise ValueError(
                    f"Expected one sample matching '{negative_sample_pattern}', found: {len(matching_negative_samples)}"
                )
            negative_sample = matching_negative_samples[0]

            # Find all samples that match the negative group pattern
            matching_group_samples = [
                sample
                for sample in sample_names
                if re.search(negative_group_pattern, sample)
            ]
            if not matching_group_samples:
                raise ValueError(
                    f"No samples found matching the group pattern '{negative_group_pattern}'"
                )

            # Add negative sample and group samples to the negative group dictionary
            negative_groups[negative_sample] = set(matching_group_samples)

    return negative_groups


def run(
    report_paths: list[str],
    results_path: str,
    negative_group_patterns: list[tuple[str, str]] | None,
):
    # Get output paths for the results
    combined_species_out, rrpm_out, tophits_out = get_output_paths(results_path)

    # Get sample names mapped to their report paths
    sample_names = get_sample_names(report_paths)

    # Mapping of sample name and taxID to their data
    sample_organism_data: dict[tuple[str, int], dict[str, Any]] = {}

    # Dictionary to store output combined species data
    combined_species_data: dict[int, dict[str, Any]] = {}

    # Dictionary for storing the number of reads for each sample
    n_reads = {}

    # For each report (and its corresponding sample name) in the folder
    for sample_name, report_path in sample_names.items():
        # Open the report file
        with open(report_path) as report_file:
            # Skip the first two lines (headers)
            report_file.readline()
            report_file.readline()
            report_reader = csv.DictReader(report_file, delimiter="\t")

            # For each row in the report file
            for row in report_reader:
                # Calculate total reads for the sample by summing root and unclassified reads
                if row["taxID"] in ["0", "1"]:
                    n_reads.setdefault(sample_name, 0)
                    n_reads[sample_name] += int(row["reads"])
                    continue

                # Skip rows that are not at the species level
                if row["rank"] != "species":
                    continue

                # The taxID of the current organism
                organism = int(row["taxID"])

                # Stored for future use in tophits output
                sample_organism_data[(sample_name, organism)] = {
                    "kmers": row["kmers"],
                    "dup": row["dup"],
                    "reads": row["reads"],
                    "cov": row["cov"],
                }

                # If the organism hasn't been recorded in the combined species data, make a new entry for it
                if combined_species_data.get(organism) is None:
                    # The entry for this organism contains the number of reads recorded for each sample
                    combined_species_data[organism] = {
                        s_name: 0 for s_name in sample_names
                    }

                    # The entry also contains the taxID, taxName and total number of reads for the organism
                    combined_species_data[organism]["taxID"] = organism
                    combined_species_data[organism]["taxName"] = row[
                        "taxName"
                    ].strip()  # damn you kraken
                    combined_species_data[organism]["Total # of Reads"] = 0

                # For the current organism and sample, add the number of reads
                combined_species_data[organism][sample_name] += int(row["reads"])
                combined_species_data[organism]["Total # of Reads"] += int(row["reads"])

    # Display the output rows and columns in alphabetical order
    sample_names = list(map(int, sample_names))
    sample_names.sort()
    sample_names = list(map(str, sample_names))
    organisms = list(combined_species_data.keys())
    organisms.sort()

    # Format data into list of dictionaries
    combined_species_list = [combined_species_data[organism] for organism in organisms]

    # Write the data to combined_species_out
    with open(combined_species_out, "w") as combined_species_file:
        writer = csv.DictWriter(
            combined_species_file,
            fieldnames=["taxID", "taxName", "Total # of Reads"] + sample_names,
        )

        writer.writeheader()

        for row in combined_species_list:
            writer.writerow({x: str(y) for x, y in row.items()})

    # Calculate RPM for each organism and sample
    rpm_data = []
    for row in combined_species_list:
        rpm_row = {
            "taxID": row["taxID"],
            "taxName": row["taxName"],
            "Total # of Reads": row["Total # of Reads"],
        }

        # For the sample columns, divide the value by num_million_reads for the sample
        for sample in sample_names:
            rpm_row[sample] = int(row[sample]) / (n_reads[sample] / 1_000_000)

        # Add new row to rpm_data
        rpm_data.append(rpm_row)

    # For each organism, calculate z scores in each sample of their rpm
    for row in rpm_data:
        # z scores for a single organism, across all samples
        z_scores = [x for x in stats.zscore([row[sample] for sample in sample_names])]

        for sample, z_score in zip(sample_names, z_scores):
            if sample_organism_data.get((sample, row["taxID"])):
                sample_organism_data[(sample, row["taxID"])]["z_score"] = z_score

    # Set up negative groups
    negative_groups = get_negative_groups(sample_names, negative_group_patterns)

    # Calculate rRPM for each organism and sample
    rrpm_data = []
    for row in rpm_data:
        rrpm_row = {
            "taxID": row["taxID"],
            "taxName": row["taxName"],
            "Total # of Reads": row["Total # of Reads"],
        }

        for sample in sample_names:
            negative_sample = ""
            for negative, group in negative_groups.items():
                if sample in group:
                    negative_sample = negative
                    break

            negative_value = 0
            for s in sample_names:
                if s == negative_sample:
                    negative_value = int(row[s])
                    break

            if negative_value == 0:
                negative_value = 1

            # Add new row to rrpm_data
            rrpm_row[sample] = int(row[sample]) / int(negative_value)

        rrpm_data.append(rrpm_row)

    # Write the data to rrpm_out
    with open(rrpm_out, "w") as rrpm_file:
        writer = csv.DictWriter(
            rrpm_file,
            fieldnames=["taxID", "taxName", "Total # of Reads"] + sample_names,
        )

        writer.writeheader()

        for row in rrpm_data:
            writer.writerow({x: str(y) for x, y in row.items()})

    # Calculate tophits for each sample
    tophits_data = []
    for sample in sample_names:
        sorted_data = sorted(rrpm_data, key=lambda x: x[sample], reverse=True)
        topfifteen = [(x["taxID"], x["taxName"], x[sample]) for x in sorted_data[0:15]]

        for i, (taxid, taxname, rrpm) in enumerate(topfifteen, start=1):
            sample_org = sample_organism_data.get((sample, taxid))

            # If number of organisms in a sample is less than 15, then the top 15 will include organisms not in the sample
            # This needs to be checked
            if sample_org is not None:
                tophits_row = {
                    "sampleName": sample,
                    "taxID": taxid,
                    "taxName": taxname,
                    "rank": i,
                    "rRPM": rrpm,
                    "kmers": sample_org["kmers"],
                    "dup": sample_org["dup"],
                    "reads": sample_org["reads"],
                    "cov": sample_org["cov"],
                    "z_score": sample_org["z_score"],
                }

                tophits_data.append(tophits_row)

    # Write the data to tophits_out
    with open(tophits_out, "w") as tophits_file:
        writer = csv.DictWriter(
            tophits_file,
            fieldnames=[
                "sampleName",
                "taxID",
                "taxName",
                "rank",
                "rRPM",
                "kmers",
                "dup",
                "reads",
                "cov",
                "z_score",
            ],
        )

        writer.writeheader()

        for row in tophits_data:
            writer.writerow({x: str(y) for x, y in row.items()})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--reports",
        required=True,
        type=str,
        nargs="+",
        help="Input KrakenUniq report files",
    )
    parser.add_argument(
        "--results",
        required=False,
        type=str,
        default="results",
        help="Directory to store the output files (default: results)",
    )
    parser.add_argument(
        "--negative-group",
        required=False,
        nargs=2,
        action="append",
        metavar=("NEGATIVE_SAMPLE_PATTERN", "NEGATIVE_GROUP_PATTERN"),
        help="Negative control group definition. Provide patterns for matching the negative sample and its group.",
    )
    args = parser.parse_args()

    run(
        report_paths=args.reports,
        results_path=args.results,
        negative_group_patterns=args.negative_group,
    )


if __name__ == "__main__":
    main()

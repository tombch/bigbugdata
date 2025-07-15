import os
import re
import csv
import logging
import argparse
from typing import Any
from pathlib import Path
from importlib.metadata import version
from scipy import stats


# Set up logging
logging.basicConfig(
    format="[%(levelname)s] %(message)s",
    level=logging.INFO,
)


def get_output_paths(results_dir: str, rank: str) -> tuple[Path, Path, Path]:
    """
    Get the output file paths for combined taxa, rrpm, and tophits.
    """

    # Create the results directory if it doesn't exist
    Path(results_dir).mkdir(parents=True, exist_ok=True)

    # Output file paths
    combined_taxa_path = Path(results_dir) / f"combined_{rank}.csv"
    rrpm_path = Path(results_dir) / f"rrpm_{rank}.csv"
    tophits_path = Path(results_dir) / f"tophits_{rank}.csv"

    return combined_taxa_path, rrpm_path, tophits_path


def get_sample_ids(report_paths: list[str]) -> dict[str, str]:
    """
    Get mapping of sample IDs to their report paths by partitioning on the last underscore.
    """

    return {
        os.path.basename(report).rpartition("_")[0]: report for report in report_paths
    }


def get_ordered_sample_ids(sample_ids: list[str]) -> list[str]:
    """
    Get sample IDs ordered by their numeric value, if possible.
    If not numeric, return them sorted as strings.
    """

    try:
        # Try to convert sample IDs to integers for sorting
        ordered_sample_ids = sorted(sample_ids, key=lambda x: int(x))
    except ValueError:
        # If conversion fails, sort them as strings
        ordered_sample_ids = sorted(sample_ids)

    return ordered_sample_ids


def get_negative_control_groups(
    sample_ids: list[str],
    group_patterns: list[tuple[str, str]] | None,
) -> dict[str, set[str]]:
    """
    Get a dictionary mapping negative control samples to their groups based on the provided patterns.
    """

    negative_groups = {}

    if group_patterns is not None:
        for negative_sample_pattern, negative_group_pattern in group_patterns:
            # Find sample that matches the negative sample pattern
            matching_negative_samples = [
                sample_id
                for sample_id in sample_ids
                if re.search(negative_sample_pattern, sample_id)
            ]
            if len(matching_negative_samples) != 1:
                raise ValueError(
                    f"Expected one sample matching '{negative_sample_pattern}', found: {len(matching_negative_samples)}"
                )
            negative_sample = matching_negative_samples[0]

            # Find all samples that match the negative group pattern
            matching_group_samples = [
                sample_id
                for sample_id in sample_ids
                if re.search(negative_group_pattern, sample_id)
            ]
            if not matching_group_samples:
                raise ValueError(
                    f"No samples found matching the group pattern '{negative_group_pattern}'"
                )

            # Add negative sample and group samples to the negative group dictionary
            negative_groups[negative_sample] = set(matching_group_samples)
            logging.info(f"Negative Control ID: {negative_sample}")
            logging.info(f"Negative Control Group: {', '.join(matching_group_samples)}")

    return negative_groups


def get_rpms(
    combined_taxa_list: list[dict[str, Any]],
    n_reads: dict[str, int],
) -> list[dict[str, Any]]:
    """
    Calculate RPM for each organism and sample based on the combined taxa data and number of reads.
    """

    rpm_data = []
    for row in combined_taxa_list:
        rpm_row = {
            "taxID": row["taxID"],
            "taxName": row["taxName"],
            "Total # of Reads": row["Total # of Reads"],
        }

        # For the sample columns, divide the value by num_million_reads for the sample
        for sample_id in n_reads:
            rpm_row[sample_id] = int(row[sample_id]) / (n_reads[sample_id] / 1_000_000)

        # Add new row to rpm_data
        rpm_data.append(rpm_row)

    return rpm_data


def get_rrpms(
    rpm_data: list[dict[str, Any]],
    sample_ids: list[str],
    negative_groups: dict[str, set[str]],
) -> list[dict[str, Any]]:
    """
    Calculate rRPM for each organism and sample based on the RPM data and negative control groups.
    """

    rrpm_data = []
    for row in rpm_data:
        rrpm_row = {
            "taxID": row["taxID"],
            "taxName": row["taxName"],
            "Total # of Reads": row["Total # of Reads"],
        }

        for sample_id in sample_ids:
            negative_control_id = ""
            for negative, group in negative_groups.items():
                if sample_id in group:
                    negative_control_id = negative
                    break

            negative_control_rpm = int(row.get(negative_control_id, 1))
            if negative_control_rpm == 0:
                negative_control_rpm = 1

            # Calculate rRPM as RPM of the sample divided by RPM of the negative control
            rrpm_row[sample_id] = int(row[sample_id]) / int(negative_control_rpm)

        rrpm_data.append(rrpm_row)

    return rrpm_data


def get_tophits(
    rrpm_data: list[dict[str, Any]],
    sample_ids: list[str],
    sample_organism_data: dict[tuple[str, int], dict[str, Any]],
    n_tophits: int,
) -> list[dict[str, Any]]:
    """
    Get the top hits for each sample based on the rRPM data.
    """

    tophits_data = []
    for sample_id in sample_ids:
        sorted_data = sorted(rrpm_data, key=lambda x: x[sample_id], reverse=True)
        tophits = [
            (x["taxID"], x["taxName"], x[sample_id]) for x in sorted_data[0:n_tophits]
        ]

        for i, (taxid, taxname, rrpm) in enumerate(tophits, start=1):
            sample_org = sample_organism_data.get((sample_id, taxid))

            # If number of organisms in a sample is less than 15, then the top n_tophits will include organisms not in the sample
            # This needs to be checked <-- TODO what does this mean?
            if sample_org is not None:
                tophits_row = {
                    "sampleName": sample_id,
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

    return tophits_data


def write_file(
    file_path: Path,
    data: list[dict[str, Any]],
    fieldnames: list[str],
):
    """
    Write the data to a CSV file at the specified path.
    """

    logging.info(f"Writing data to {file_path}")
    with open(file_path, "w") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        for row in data:
            writer.writerow({x: str(y) for x, y in row.items()})


def run(
    report_paths: list[str],
    results_path: str,
    rank: str,
    n_tophits: int,
    group_patterns: list[tuple[str, str]] | None,
):
    # Output paths for the results
    combined_taxa_path, rrpm_path, tophits_path = get_output_paths(results_path, rank)

    # Mapping of sample IDs to their report paths
    sample_ids = get_sample_ids(report_paths)

    # Mapping of sample ID and taxID to their data
    sample_organism_data: dict[tuple[str, int], dict[str, Any]] = {}

    # Dictionary to store combined taxa data
    combined_taxa_data: dict[int, dict[str, Any]] = {}

    # Mapping of sample ID to total number of reads
    n_reads = {}

    # For each sample ID (and its corresponding report)
    for sample_id, report_path in sample_ids.items():
        # Open the report file
        with open(report_path) as report_file:
            # Skip the first two lines (headers)
            report_file.readline()
            report_file.readline()
            report_reader = csv.DictReader(report_file, delimiter="\t")

            # For each row in the report file
            for row in report_reader:
                # Calculate total reads for the sample by adding root and unclassified read counts
                if row["taxID"] in ["0", "1"]:
                    n_reads.setdefault(sample_id, 0)
                    n_reads[sample_id] += int(row["reads"])
                    continue

                # Skip all other rows that are not at the specified rank
                if row["rank"] != rank:
                    continue

                # The taxID of the current organism
                organism = int(row["taxID"])

                # Stored for future use in tophits output
                sample_organism_data[(sample_id, organism)] = {
                    "kmers": row["kmers"],
                    "dup": row["dup"],
                    "reads": row["reads"],
                    "cov": row["cov"],
                }

                # If the organism hasn't been recorded in the combined taxa data, make a new entry for it
                if combined_taxa_data.get(organism) is None:
                    # The entry for this organism contains the number of reads recorded for each sample
                    combined_taxa_data[organism] = {
                        sample_id: 0 for sample_id in sample_ids
                    }

                    # The entry also contains the taxID, taxName and total number of reads for the organism
                    combined_taxa_data[organism]["taxID"] = organism
                    combined_taxa_data[organism]["taxName"] = row[
                        "taxName"
                    ].strip()  # damn you kraken
                    combined_taxa_data[organism]["Total # of Reads"] = 0

                # For the current organism and sample, add the number of reads
                combined_taxa_data[organism][sample_id] += int(row["reads"])
                combined_taxa_data[organism]["Total # of Reads"] += int(row["reads"])

    # Display the output rows and columns in alphabetical order
    sample_ids = get_ordered_sample_ids(list(sample_ids))

    # Format data into list of dictionaries
    combined_taxa_list = [
        combined_taxa_data[organism] for organism in sorted(list(combined_taxa_data))
    ]

    # Write the combined_taxa file
    write_file(
        file_path=combined_taxa_path,
        data=combined_taxa_list,
        fieldnames=["taxID", "taxName", "Total # of Reads"] + sample_ids,
    )

    # Calculate RPM for each organism and sample
    rpm_data = get_rpms(combined_taxa_list, n_reads)

    # For each organism, calculate z scores in each sample of their rpm
    for row in rpm_data:
        # z scores for a single organism, across all samples
        z_scores = [
            x for x in stats.zscore([row[sample_id] for sample_id in sample_ids])
        ]

        for sample_id, z_score in zip(sample_ids, z_scores):
            if sample_organism_data.get((sample_id, row["taxID"])):
                sample_organism_data[(sample_id, row["taxID"])]["z_score"] = z_score

    # Set up negative groups
    negative_groups = get_negative_control_groups(sample_ids, group_patterns)

    # Calculate rRPM for each organism and sample
    rrpm_data = get_rrpms(rpm_data, sample_ids, negative_groups)

    # Write the rrpm file
    write_file(
        file_path=rrpm_path,
        data=rrpm_data,
        fieldnames=["taxID", "taxName", "Total # of Reads"] + sample_ids,
    )

    # Calculate tophits for each sample
    tophits_data = get_tophits(rrpm_data, sample_ids, sample_organism_data, n_tophits)

    # Write the tophits file
    write_file(
        file_path=tophits_path,
        data=tophits_data,
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-r",
        "--reports",
        required=True,
        type=str,
        nargs="+",
        help="Input KrakenUniq report files",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=False,
        type=str,
        default="results",
        help="Directory to store the output files (default: results)",
    )
    parser.add_argument(
        "-n",
        "--nc-group",
        required=False,
        nargs=2,
        action="append",
        metavar=("CONTROL", "GROUP"),
        help="Provide REGEX patterns to match a negative control and its group of samples",
    )
    parser.add_argument(
        "-R",
        "--rank",
        required=False,
        type=str,
        default="species",
        help="Taxonomic rank to filter the reports by (default: species)",
    )
    parser.add_argument(
        "-t",
        "--tophits",
        required=False,
        type=int,
        default=15,
        help="Number of top hits to include in the tophits output (default: 15)",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {version('bigbugdata')}",
        help="Output the version of bigbugdata",
    )

    args = parser.parse_args()
    logging.info(f"bigbugdata v{version('bigbugdata')}")

    run(
        report_paths=args.reports,
        results_path=args.output,
        rank=args.rank,
        n_tophits=args.tophits,
        group_patterns=args.nc_group,
    )


if __name__ == "__main__":
    main()

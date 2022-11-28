from scipy import stats
import argparse
import glob
import csv
import os


def run(
    species_level_reports,
    dna_totalreads,
    rna_totalreads,
    combined_species_out,
    rrpm_out,
    tophits_out,
):
    # List of file paths for each input kraken tsv
    species_level_paths = glob.glob(f"{species_level_reports}/*.tsv")

    # Sample names from each tsv
    sample_names = [
        os.path.basename(tsv_path).split("_")[0] for tsv_path in species_level_paths
    ]

    # Data needed for final tophits output
    sample_organism_data = {}

    # Dictionary to store output combined species data
    combined_species_data = {}

    # For each tsv (and its corresponding sample name) in the folder
    for tsv_path, sample_name in zip(species_level_paths, sample_names):

        # Open the tsv file
        with open(tsv_path) as tsv_file:
            tsv_reader = csv.DictReader(tsv_file, delimiter="\t")

            # For each row in the tsv file
            for row in tsv_reader:
                # The name of the current organism
                organism = row["taxName"]

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

                    # The entry also contains the name of the organism
                    combined_species_data[organism]["taxName"] = organism

                    # The entry also contains the total number of reads for the organism
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

    # Reformat data into list of dictionaries
    combined_species_data = [combined_species_data[organism] for organism in organisms]

    # Write the data to combined_species_out
    with open(combined_species_out, "w") as combined_species_file:
        writer = csv.DictWriter(
            combined_species_file,
            fieldnames=["taxName", "Total # of Reads"] + sample_names,
        )

        writer.writeheader()

        for row in combined_species_data:
            writer.writerow(row)

    # Dictionary for storing the num million reads for each sample
    num_million_reads = {}

    # Open DNA tsv and add the sample number + num million reads for each sample
    with open(dna_totalreads) as dna_file:
        dna_data = csv.reader(dna_file, delimiter="\t")
        for row in dna_data:
            sample = row[0]
            num = row[3]
            num_million_reads[sample.split("_")[0]] = float(num)

    # Open RNA tsv and add the sample number + num million reads for each sample
    with open(rna_totalreads) as rna_file:
        rna_data = csv.reader(rna_file, delimiter="\t")
        for row in rna_data:
            sample = row[0]
            num = row[3]
            num_million_reads[sample.split("_")[0]] = float(num)

    # Calculate RPM for each organism and sample
    rpm_data = []
    for row in combined_species_data:
        rpm_row = {
            "taxName": row["taxName"],
            "Total # of Reads": row["Total # of Reads"],
        }

        # For the sample columns, divide the value by num_million_reads for the sample
        for sample in sample_names:
            rpm_row[sample] = int(row[sample]) / num_million_reads[sample]

        # Add new row to rpm_data
        rpm_data.append(rpm_row)

    for row in rpm_data:
        z_scores = stats.zscore([row[str(x)] for x in range(1, 96)])
        for i, value in enumerate(z_scores, start=1):
            if sample_organism_data.get((str(i), row["taxName"])):
                sample_organism_data[(str(i), row["taxName"])]["z_score"] = value

    # Set up negative groups
    negative_group = {
        "6": {str(x) for x in range(1, 7)},
        "18": {str(x) for x in range(7, 19)},
        "30": {str(x) for x in range(19, 31)},
        "42": {str(x) for x in range(31, 43)},
        "48": {str(x) for x in range(43, 49)},
        "54": {str(x) for x in range(49, 55)},
        "66": {str(x) for x in range(55, 67)},
        "78": {str(x) for x in range(67, 79)},
        "90": {str(x) for x in range(79, 91)},
        "95": {str(x) for x in range(91, 96)},
    }

    # Calculate rRPM for each organism and sample
    rrpm_data = []
    for row in rpm_data:
        rrpm_row = {
            "taxName": row["taxName"],
            "Total # of Reads": row["Total # of Reads"],
        }

        for sample in sample_names:
            negative_sample = ""
            for negative, group in negative_group.items():
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
            fieldnames=["taxName", "Total # of Reads"] + sample_names,
        )

        writer.writeheader()

        for row in rrpm_data:
            writer.writerow(row)

    # Calculate tophits for each sample
    tophits_data = []
    for sample in sample_names:
        sorted_data = sorted(rrpm_data, key=lambda x: x[sample], reverse=True)
        topfifteen = [(x["taxName"], x[sample]) for x in sorted_data[0:15]]

        for i, (taxname, rrpm) in enumerate(topfifteen, start=1):
            sample_org = sample_organism_data.get((sample, taxname))

            # If number of organisms in a sample is less than 15, then the top 15 will include organisms not in the sample
            # This needs to be checked
            if sample_org is not None:
                tophits_row = {
                    "sampleName": sample,
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
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--species-level-reports",
        required=True,
        help="Input directory of species level report .tsv files",
    )
    parser.add_argument(
        "--dna-totalreads", required=True, help="Input .tsv file of DNA total reads"
    )
    parser.add_argument(
        "--rna-totalreads", required=True, help="Input .tsv file of RNA total reads"
    )

    parser.add_argument(
        "--combined-species-out",
        required=True,
        help="Output .csv file for the combined species data",
    )
    parser.add_argument(
        "--rrpm-out", required=True, help="Output .csv file for the rrpm data"
    )
    parser.add_argument(
        "--tophits-out",
        required=True,
        help="Output .csv file for the tophits data",
    )

    args = parser.parse_args()

    run(
        species_level_reports=args.species_level_reports,
        dna_totalreads=args.dna_totalreads,
        rna_totalreads=args.rna_totalreads,
        combined_species_out=args.combined_species_out,
        rrpm_out=args.rrpm_out,
        tophits_out=args.tophits_out,
    )


if __name__ == "__main__":
    main()

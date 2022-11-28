import argparse
import csv
import glob


def run(
    kraken_reports,
    dna_totalreads,
    rna_totalreads,
    combined_kraken_out,
    rrpm_out,
    tophits_out,
):
    # List of file paths for each input kraken tsv
    kraken_tsv_paths = glob.glob(f"{kraken_reports}/*.tsv")

    # Sample names for each tsv
    sample_names = [
        tsv_path.split("/")[-1].split("_")[0] for tsv_path in kraken_tsv_paths
    ]

    # Data needed for final tophits output
    sample_organism_kmer_dup_data = {}

    # Dictionary to store output combined kraken data
    combined_kraken_data = {}

    # For each csv (and its corresponding sample name) in the folder
    for tsv_path, sample_name in zip(kraken_tsv_paths, sample_names):

        # Open the csv file
        with open(tsv_path) as tsv_file:
            tsv_reader = csv.DictReader(tsv_file, delimiter="\t")

            # For each row in the tsv file
            for row in tsv_reader:
                # The name of the current organism
                organism = row["taxName"]

                # Stored for future use in tophits output
                sample_organism_kmer_dup_data[(sample_name, organism)] = (
                    row["kmers"],
                    row["dup"],
                )

                # If the organism hasn't been recorded in the combined kraken data, make a new entry for it
                if combined_kraken_data.get(organism) is None:

                    # The entry for this organism contains the number of reads recorded for each sample
                    combined_kraken_data[organism] = {
                        s_name: 0 for s_name in sample_names
                    }

                    # The entry also contains the name of the organism
                    combined_kraken_data[organism]["taxName"] = organism

                    # The entry also contains the total number of reads for the organism
                    combined_kraken_data[organism]["Total # of Reads"] = 0

                # For the current organism and sample, add the number of reads
                combined_kraken_data[organism][sample_name] += int(row["reads"])
                combined_kraken_data[organism]["Total # of Reads"] += int(row["reads"])

    # Display the output rows and columns in alphabetical order
    sample_names.sort()
    organisms = list(combined_kraken_data.keys())
    organisms.sort()

    # Write the data to combined_kraken_out
    with open(combined_kraken_out, "w") as combined_kraken_file:
        writer = csv.DictWriter(
            combined_kraken_file,
            fieldnames=["taxName", "Total # of Reads"] + sample_names,
        )
        writer.writeheader()
        for organism in organisms:
            writer.writerow(combined_kraken_data[organism])

    # Re-open csv and add each row
    combined_kraken_rows = []
    with open(combined_kraken_out, "r") as combined_kraken_file:
        data = csv.reader(combined_kraken_file)
        combined_kraken_header = next(data)
        for row in data:
            combined_kraken_rows.append(row)

    # Dictionary for storing the num million reads for each sample
    num_million_reads = {}

    # Open DNA tsv and add the sample number + num million reads for each sample
    with open(dna_totalreads) as dna_file:
        data = csv.reader(dna_file, delimiter="\t")
        for row in data:
            if row:
                sample = row[0]
                num = row[3]
                num_million_reads[sample.split("_")[0]] = float(num)

    # open RNA tsv and add the sample number + num million reads for each sample
    with open(rna_totalreads) as rna_file:
        data = csv.reader(rna_file, delimiter="\t")
        for row in data:
            if row:
                sample = row[0]
                num = row[3]
                num_million_reads[sample.split("_")[0]] = float(num)

    # Store the order of sample numbers from the csv in a list
    samples = [column.split("_")[0] for column in combined_kraken_header[2:]]

    rows_rpm = []
    for row in combined_kraken_rows:
        new_row = []

        # Add the first two columns of data, unchanged
        for value in row[0:2]:
            new_row.append(value)

        # For the sample columns, divide the value by num_million_reads for the sample
        for sample, value in zip(samples, row[2:]):
            rpm = int(value) / num_million_reads[sample]
            new_row.append(rpm)

        # Add new row to rows_rpm
        rows_rpm.append(new_row)

    # Set up dictionary of sets
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

    rRPM = []
    for row in rows_rpm:
        new_row = []

        # Add the first two columns of data, unchanged
        for value in row[0:2]:
            new_row.append(value)

        for sample, value in zip(samples, row[2:]):
            negative_sample = ""
            for negative, group in negative_group.items():
                if sample in group:
                    negative_sample = negative
                    break

            negative_value = 0
            for s, v in zip(samples, row[2:]):
                if s == negative_sample:
                    negative_value = int(v)
                    break

            if negative_value == 0:
                negative_value = 1

            rrpm = int(value) / int(negative_value)
            new_row.append(rrpm)

        rRPM.append(new_row)

    # Write the output to a file
    with open(rrpm_out, "w") as fh_out:
        writer = csv.writer(fh_out)
        writer.writerow(combined_kraken_header)
        for row in rRPM:
            writer.writerow(row)

    with open(tophits_out, "w") as tophits_out_file:
        writer = csv.writer(tophits_out_file)
        writer.writerow(["taxName", "rRPM", "kmers", "dup"])

        # top15
        for n in range(len(samples)):
            sorted_data = sorted(rRPM[1:], key=lambda x: x[n + 2], reverse=True)
            topfifteen = [(x[0], x[n + 2]) for x in sorted_data[0:16]]
            writer.writerow(samples[n])
            for taxname, rrpm in topfifteen:
                kmers, dup = sample_organism_kmer_dup_data[(samples[n], taxname)]
                writer.writerow([taxname, rrpm, kmers, dup])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--kraken-reports", required=True, help="Input directory of kraken reports"
    )
    parser.add_argument(
        "--dna-totalreads", required=True, help="Input tsv of DNA total reads"
    )
    parser.add_argument(
        "--rna-totalreads", required=True, help="Input tsv of RNA total reads"
    )

    parser.add_argument(
        "--combined-kraken-out",
        required=True,
        help="Output filename for the combined kraken data",
    )
    parser.add_argument(
        "--rrpm-out", required=True, help="Output file name for the rrpm data"
    )
    parser.add_argument(
        "--tophits-out", required=True, help="Output file name for the tophits data"
    )

    args = parser.parse_args()

    run(
        kraken_reports=args.kraken_reports,
        dna_totalreads=args.dna_totalreads,
        rna_totalreads=args.rna_totalreads,
        combined_kraken_out=args.combined_kraken_out,
        rrpm_out=args.rrpm_out,
        tophits_out=args.tophits_out,
    )


if __name__ == "__main__":
    main()

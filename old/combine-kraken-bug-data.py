import csv
import sys
import glob
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("folder")
    args = parser.parse_args()

    # List of file paths for each input csv
    csv_paths = glob.glob(f"{args.folder}/*.tsv")
    
    # Sample names for each csv
    sample_names = [csv_path.split("/")[-1].split("_")[0] for csv_path in csv_paths]
    
    # Dictionary to store output data
    output_data = {}

    # For each csv (and its corresponding sample name) in the folder
    for csv_path, sample_name in zip(csv_paths, sample_names):

        # Open the csv file
        with open(csv_path) as csv_file:
            csv_reader = csv.DictReader(csv_file, delimiter="\t")
            
            # For each row in the csv file
            for row in csv_reader:
                # The name of the current organism
                organism = row["taxName"]

                # If the organism hasn't been recorded in the output data, make a new entry for it
                if output_data.get(organism) is None:

                    # The entry for this organism contains the number of reads recorded for each sample
                    output_data[organism] = {s_name : 0 for s_name in sample_names}

                    # The entry also contains the name of the organism
                    output_data[organism]["taxName"] = organism

                    # The entry also contains the total number of reads for the organism
                    output_data[organism]["Total # of Reads"] = 0

                #For the current organism and sample, add the number of reads w/ children
                output_data[organism][sample_name] += int(row["reads"])
                output_data[organism]["Total # of Reads"] += int(row["reads"])
    
    # Display the output rows and columns in alphabetical order
    sample_names.sort()
    organisms = list(output_data.keys())
    organisms.sort()

    # Write the data to stdout
    writer = csv.DictWriter(sys.stdout, fieldnames=["taxName", "Total # of Reads"] + sample_names)
    writer.writeheader()
    for organism in organisms:
        writer.writerow(output_data[organism])
                    

if __name__ == '__main__':
    main()
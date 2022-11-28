import csv

# Open csv and add each row to rows
rows = []
with open('meningoencephalitis/Results/Kraken_report/combine-kraken-bug-data.csv', 'r') as fh_in:
    data = csv.reader(fh_in)
    header = next(data)
    for row in data:
        rows.append(row)

#print number of rows
print(len(rows))       

#For every row, check it meets the filter criteria and if it does add to filt
#filt = []
#for row in rows:
 #    if any(i > 5 for i in list(map(int, row[2:]))):
  #      filt.append(row)

#Write the filtered output to a file
#with open('meningoencephalitis/Results/DNA_filtering/filtered_combined_DNA_bugdata_19-08-22.csv', 'w') as fh_out:
 #    writer = csv.writer(fh_out)
  #   writer.writerow(header)
   #  for row in filt:
    #     writer.writerow(row)


# Dictionary for storing the num million reads for each sample
num_million_reads = {}

# Open tsv and add the sample number + num million reads for each sample
with open('meningoencephalitis/DNA_totalreads.tsv') as dna_file:
    data = csv.reader(dna_file, delimiter='\t')
    for row in data:
        if row:
            sample = row[0]
            num = row[3]
            num_million_reads[sample.split('_')[0]] = float(num)

#open RNA tsv and add the 
with open('meningoencephalitis/RNA_totalreads_sort.tsv') as dna_file:
    data = csv.reader(dna_file, delimiter='\t')
    for row in data:
        if row:
            sample = row[0]
            num = row[3]
            num_million_reads[sample.split('_')[0]] = float(num)

# Store the order of sample numbers from the csv in a list
samples = [column.split('_')[0] for column in header[2:]]

rows_rpm = []
for row in rows:
    new_row = []

    # Add the first two columns of data, unchanged
    for value in row[0:2]:
        new_row.append(value)

    # For the sample columns, divide the value by num_million_reads for the sample
    for sample, value in zip(samples, row[2:]):
        rpm = int(value) / num_million_reads[sample]
        new_row.append(rpm)
    
    # Add new row to filt_with_rpm
    rows_rpm.append(new_row)

# Write the output to a new file
with open('meningoencephalitis/Results/Kraken_report/RPM_Krakenout_23-09-22.csv', 'w') as fh_out:
    writer = csv.writer(fh_out)
    writer.writerow(header)
    for row in rows_rpm:
        writer.writerow(row)

# print
#print(len(rows_rpm))

# Filter RPM less than 10
#filt_RPM = []
#for row in filt_with_rpm:
 #   if any(i > 10 for i in list(map(int, row[2:]))):
  #      filt_RPM.append(row)
 
# print
#print(len(filt_RPM))

# Write the filtered output to a file
#with open('meningoencephalitis/Results/DNA_filtering/filtered_RPM_19-08-22.csv', 'w')as fh_out:
 #   writer = csv.writer(fh_out)
  #  writer.writerow(header)
   # for row in filt_RPM:
    #    writer.writerow(row)

# Set up dictionary of sets
negative_group = {
    "6" : {str(x) for x in range(1, 7)},
    "18" : {str(x) for x in range(7, 19)},
    "30" : {str(x) for x in range(19, 31)},
    "42" : {str(x) for x in range(31, 43)},
    "48" : {str(x) for x in range(43, 49)},
    "54" : {str(x) for x in range(49, 55)},
    "66" : {str(x) for x in range(55, 67)},
    "78" : {str(x) for x in range(67, 79)},
    "90" : {str(x) for x in range(79, 91)},
    "95" : {str(x) for x in range(91, 96)}
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
with open('meningoencephalitis/Results/Kraken_report/rRPM_Krakenout_23-09-22.csv', 'w')as fh_out:
    writer = csv.writer(fh_out)
    writer.writerow(header)
    for row in rRPM:
        writer.writerow(row)

#Filter rRPM <10
#filt_rRPM = []
#for row in rRPM:
 #   if any(i > 10 for i in list(map(int, row[2:]))):
  #      filt_rRPM.append(row)

#print length
#print(len(filt_rRPM))

# Write the filtered output to a file
#with open('meningoencephalitis/Results/Kraken_report/filt_rRPM_Krakenout_23-09-22.csv', 'w')as fh_out:
 #   writer = csv.writer(fh_out)
  #  writer.writerow(header)
   # for row in filt_rRPM:
    #    writer.writerow(row)

#top15
for n in range(len(samples)):
    sorted_data = sorted(rRPM[1:], key=lambda x:x[n+2], reverse=True)
    topfifteen = [(x[0],x[n+2]) for x in sorted_data[0:16]]
    print(samples[n])
    for top in topfifteen:
        print(top)


"""
df = pd.read_csv('meningoencephalitis/Results/combined_reads_07-07-22.csv', index_col="Organism Name")
#filter all rows in which the total numnber of reads for that organism is less than or equal to 3
df_filtered = df[df['Total # of Reads (w/ Children)'] > 4]
#write filtered filt to csv
df_filtered.to_csv('meningoencephalitis/Results/filtered_combined_reads_07-07-22.csv',index=True)
"""
# GTF formatting and filtering utility for microbial reference genomes.
# This script is part of the genome preprocessing pipeline (MGBC/MGnify/RefSeq)
# Format gtf file
def convert_gtf(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            id_index = line.find('ID ')
            if id_index != -1:
                # Replace space to tab
                before_id = line[:id_index].replace(' ', '\t')
                after_id = line[id_index:]  # Keep the original content after "ID "
                outfile.write(before_id + after_id)
            else:
                outfile.write(line)

def filter_records(input_file, exclusion_file, included_output, excluded_output):
    with open(exclusion_file, 'r') as ef:
        exclusion_ids = set(line.strip() for line in ef)

    # Set output
    with open(input_file, 'r') as infile, \
         open(included_output, 'w') as included_file, \
         open(excluded_output, 'w') as excluded_file:
        
        for line in infile:
            # Write to excluded file if line contains ncRNA
            if "ncRNA" in line:
                excluded_file.write(line)
                continue
                
            if "ID " in line:
                id_start = line.index("ID ") + 3
                id_end = line.index(";", id_start) if ";" in line[id_start:] else len(line)
                record_id = line[id_start:id_end].strip('"')

                # Write to corresponding output file
                if record_id in exclusion_ids:
                    excluded_file.write(line)
                else:
                    included_file.write(line)


convert_gtf('M1_CDS_ncRNA.gtf', 'M1_genome_all.fix.gtf')

filter_records(
    input_file="M1_genome_all.fix.gtf",
    exclusion_file="null.txt",  # File containing the IDs to exclude, empty file means no exclusion
    included_output="M1_genome_all.cds.gtf",  # Output file for records in the set
    excluded_output="M1_genome_all_ncRNA.gtf"   # Output file for records not in the set
)

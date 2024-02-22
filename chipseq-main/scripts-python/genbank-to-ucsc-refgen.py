#!/usr/bin/python3
"""Usage: refseq-to-ucsc-refgen.py [1: refgen.fasta] [2: report_url] [3: output/path]"""
"""Description: looks through your target reference genome (refgen.fasta) and replaces GenBank annotations with
   UCSC annotations based on the URL of the assembly report as found on NCBI (or whatever you're using)."""

### define url scraper to gain the assembly report_path
def parse_assembly_report(report_url):
    ## Assertions
    assert report_url.endswith('.txt'), f"The url {report_url} does not end with '.txt'. Exiting..."

    ## Init assembly report
    assembly_data = []
    with urllib.request.urlopen(report_url) as report:
        for line in report:
            line = line.decode().strip()
            if line.startswith('#'):
                if line.startswith('# Sequence-Name'):
                    header = line[1:].strip().split('\t')
                continue
            fields = line.split('\t')
            assembly_data.append(fields)

    ## Convert to DataFrame
    matcher = pd.DataFrame(assembly_data, columns=header)

    # Keep only the RefSeq-Accn and UCSC-style-name columns
    matcher = matcher[['GenBank-Accn', 'UCSC-style-name']]

    return matcher


### define read_genome to read the input reference genome in a memory-efficient way (more memory efficient than loading the whole thing into memory)
def read_genome(file_path):
    genome = {}
    with open(file_path, 'r', encoding='ISO-8859-1') as infile:
        seq_id = None
        seq = []
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if seq_id is not None:
                    genome[seq_id] = ''.join(seq)
                    seq = []
                seq_id = line[1:]
            else:
                seq.append(line)
        genome[seq_id] = ''.join(seq)
    return genome


### define the reference genome substitutor
def genome_sub(refgen_path, output_path, matcher):
        # init genome_counter
        genome_counter = 0
        # Read in the reference genome
        genome = read_genome(refgen_path)
        # Run substitution loop
        with open(output_path, 'w') as outfile:
            if refgen_path.endswith('.gz'):
                infile = gzip.open(refgen_path, 'rt')
            else:
                infile = open(refgen_path, 'r')
            seq_id = None
            seq = []
            for line in infile:
                line = line.strip()
                if line.startswith('>'):
                    if seq_id is not None:
                        # loop through each row of matcher and replace all occurrences of GenBank-Accn with UCSC-style-name
                        new_seq_id = seq_id
                        for i, row in matcher.iterrows():
                            new_seq_id = new_seq_id.replace(row['GenBank-Accn'], row['UCSC-style-name'])
                        outfile.write('>' + new_seq_id + '\n')
                        outfile.write(''.join(seq) + '\n')
                        seq = []
                        genome_counter += 1
                    seq_id = line[1:]
                else:
                    seq.append(line)
            # write the last sequence to output
            new_seq_id = seq_id
            for i, row in matcher.iterrows():
                new_seq_id = new_seq_id.replace(row['GenBank-Accn'], row['UCSC-style-name'])
            outfile.write('>' + new_seq_id + '\n')
            outfile.write(''.join(seq) + '\n')
        # make sure the opened infile is closed
        infile.close()
        # print number of modifications
        print("Reference genome modification finished! ", genome_counter, " lines have been replaced.")


### define the annotation substitutor
def gff_sub(gff_path, output_path, matcher):
    with gzip.open(gff_path, 'rt') as infile:
        with open(output_path.replace('.fasta', '.gff'), 'w') as outfile:
            gff_counter = 0
            for line in infile:
                if line.startswith('#'):
                    outfile.write(line)
                else:
                    fields = line.strip().split('\t')
                    seq_id = fields[0]
                    start = fields[3]
                    end = fields[4]
                    strand = fields[6]
                    attributes = fields[8]
                    # loop through each row of matcher and replace all occurrences of GenBank-Accn with UCSC-style-name
                    for i, row in matcher.iterrows():
                        seq_id = seq_id.replace(row['GenBank-Accn'], row['UCSC-style-name'])
                    # write the modified line to output
                    outfile.write('\t'.join([seq_id, fields[1], fields[2], start, end, fields[5], strand, fields[7], attributes]) + '\n')
                    gff_counter += 1
    # print number of modifications
    print("GFF modification finished! ", gff_counter, " lines have been replaced.")


### define main
def main(refgen_path, gff_path, report_url, output_path):
    ## Load the assembly report
    try:
        matcher = parse_assembly_report(report_url)
    except:
        print("The inputted assembly report URL cannot be used by this script. This a problem with the script, as the URL seems fine.")
        sys.exit(1)

    # init replacement counters
    genome_counter = 0

    # Loop through each line of the genome and change each line matching matcher['RefSeq-Accn'] with the string in the corresponding row in matcher['UCSC-style-name']
    if gff_path == None:
        # run reference genome substitution loop
        genome_sub(refgen_path, output_path, matcher)
    else:
        # run genome annotation substitution loop
        gff_sub(gff_path, output_path, matcher)

### Main
if __name__ == '__main__':
    ## Imports
    import pandas as pd
    import sys
    import urllib.request
    import urllib.error
    import gzip

    ## Syntax check
    if len(sys.argv) != 4: # includes the script name: +1
        print("Usage: genbank-to-ucsc-refgen.py <refgen.fasta> <report_url> <output/path>")
        sys.exit(1)

    ## Get command-line arguments and check format
    if sys.argv[1].endswith('.gff.gz'):
        gff_path = sys.argv[1]
        refgen_path = None
    else:
        assert refgen_path.endswith('.fasta') or refgen_path.endswith('.fna.gz'), f"{refgen_path} is not in .fasta/.fna.gz/.gff.gz format. Either change it or this script."
        refgen_path = sys.argv[1]
        gff_path = None
    report_url = sys.argv[2]
    output_path = sys.argv[3]

    ## url check
    try:
        response = urllib.request.urlopen(report_url)
    except urllib.error.URLError as e:
        print("Error opening URL:", e.reason)

    ## run main
    main(refgen_path, gff_path, report_url, output_path)

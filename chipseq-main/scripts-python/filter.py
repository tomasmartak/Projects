#!/usr/bin/python3
"""Usage: filter.py [1: toExcise.txt] [2: input_path.fasta] [3: output_path.fasta]"""
"""Description: removes chromosomes specified in 'toExcise.txt' from 'input_path.fasta'. Make sure that the naming conventions are consistent in the two. """

"""
Additional comments: toExcise.txt: a text file that contains a list of different chromosomes you want to removed
                                   (using syntax used in 'input_path.fasta'), each name in a new line, each preceded by '>'.
                     input_path.fasta: the .fasta file you want to modify.
"""

### Defs
def list_ids(ids):
    """
    Convert the provided toExcise.txt file into a
    set of identifiers, line by line, starting with ">"
    """

    # read the first file given and generate a set (faster iteration respect lists
    identifiers = set([]) # init identifiers

    with open(ids, 'r') as fi: # loads the first script argument
        for line in fi:
            line = line.strip()
            identifiers.add(str(line).replace(">", ""))
    # return set
    return identifiers

def filter(ids, infile, outfile):
    """
    Writes a file containing only the sequences with identifiers NOT
    present in the 'identifiers' set
    """
    # assertions
    assert ids.endswith('.txt'), f"{ids} is not in .txt format"
    assert infile.endswith('.fasta'), f"{infile} is not in .fasta format"
    assert outfile.endswith('.fasta'), f"{outfile} is not in .fasta format"

    # filter out the identifiers from the
    identifiers = list_ids(ids)

    with open(infile) as original_fasta, open(outfile, 'w') as corrected_fasta: # loads the 2nd and 3rd script arguments
        records = SeqIO.parse(original_fasta, 'fasta')
        for record in records:
            print(record.id)
            if record.id not in identifiers:
                SeqIO.write(record, corrected_fasta, 'fasta')

### Code
if __name__ == '__main__':
    # Imports
    import sys
    from Bio import SeqIO

    # syntax check
    if len(sys.argv) != 3:
        print("Usage: python3 filter.py <toExcise.txt> <input_path.fasta> <output_path.fasta>")
        sys.exit(1)

    ## Code
    ids, infile, outfile = sys.argv[1:]
    filter(ids, infile, outfile)
    print("Finished run with diff.py successfully!")

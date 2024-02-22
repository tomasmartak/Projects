#!/usr/bin/python3
"""Usage: deff.py [1: largerFile.fasta] [2: smallerFile.fasta] [3: output]"""
"""Description: compares two fasta files and outputs the product of the first (hopefully larger) one minus the second (hopefully smaller) one."""

### define a file reader
def read_fasta(file):
    header = None
    sequence = []
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(sequence)
                header = line
                sequence = []
            else:
                sequence.append(line)
    if header is not None:
        yield header, ''.join(sequence)

### define a file writer
def write_fasta(data, file):
    with open(file, 'w') as f:
        for header, sequence in data:
            f.write(header + '\n')
            f.write(sequence + '\n')

### define a function that produces the difference between file 1 and file 2
def main(file_1, file_2, outfile):
    assert file_1.endswith('.fasta'), f"{file_1} is not in .fasta format"
    assert file_2.endswith('.fasta'), f"{file_2} is not in .fasta format"
    assert outfile.endswith('.fasta'), f"{outfile} is not in .fasta format"

    a = dict(read_fasta(file_1))
    b = set(read_fasta(file_2))

    # Get the part of file_1 that isn't in file_2
    c = [(header, sequence) for header, sequence in a.items() if (header, sequence) not in b]

    # Write the result to a fasta file
    write_fasta(c, outfile)

### Code
if __name__ == '__main__':
    ## Imports
    import sys

    ## syntax check
    if len(sys.argv) != 4:
        print("Usage: python3 diff.py file_1 file_2 outfile")
        sys.exit(1)

    file_1, file_2, outfile = sys.argv[1:]
    main(file_1, file_2, outfile)
    print("Finished run with diff.py successfully!")

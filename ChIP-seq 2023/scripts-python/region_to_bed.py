#!/usr/bin/python3
"""Usage: region_to_bed.py [1: chromosome] [2: start] [3: end] [4: outfile]"""
"""Description: Converts your input into a valid bed format. Take care with chromosome naming conventions."""

### Defs
## make a pandas dataframe
def makeframe(chrom, start, end):
    # assertions
    assert (' ' in chrom) == False, f"{chrom} contains a space. Please remove it"
    assert start.isnumeric() == True, f"{start} isn't a numberic. Please insert a valid numeric region without commas."
    assert end.isnumeric() == True, f"{end} isn't a numberic. Please insert a valid numeric region without commas."

    # setup pandas dataframe
    tempframe = {'Chromosome': [str(chrom)], 'Start': [int(start)], 'End': [int(end)]}
    tempframe = pd.DataFrame(tempframe)
    df = pd.concat([df, tempframe], ignore_index=True)

## main
def main(a, b, c, d):
    # assertions
    assert (' ' in a) == False, f"{a} contains a space. Please remove it"
    assert b.isnumeric() == True, f"{b} isn't a number. Please insert a valid integer without commas."
    assert c.isnumeric() == True, f"{c} isn't a number. Please insert a valid integer without commas."

    # setup pandas dataframe
    df = {'Chromosome': [str(a)], 'Start': [int(b)], 'End': [int(c)]}
    df = pd.DataFrame(df)

    # loop makeframe --> commented out; this cannot run in RMarkdown, as you cannot input into it
    check = 'check' # setup dummy variable
    while check != 'c':
        check = input("Input, separated by a single space each, [1] Chromosome [2] Start [3] End. Alternatively, pres 'c' to continue.\n")
        try:
            chrom = check.split()[0]
            start = check.split()[1]
            end = check.split()[2]
            makeframe(chrom, start, end)
        except:
            break

    # convert df to bed
    bf = pybed.BedFrame.from_frame(meta=[], data=df)
    try:
        bf.to_file(d)
        print("Finished run with region_to_bed.py successfully!")
    except:
        print("There is something wrong with your filepath. However, here is what your table would look like:")
        print(bf)

### Code
if __name__ == '__main__':
    # imports
    import pandas as pd
    from fuc import pybed
    import sys
    import os

    # syntax check
    if len(sys.argv) != 5:
        print("Incorrect number of arguments. Usage: python3 region_to_bed.py chrom start end outfile.")
        sys.exit(1)

    # definitions
    chrom, start, end, outfile = sys.argv[1:5]

    # code
    main(chrom, start, end, outfile)

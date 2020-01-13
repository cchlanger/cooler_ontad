import pandas as pd
import argparse


p = argparse.ArgumentParser(description='Convert to bedpe',
                            prefix_chars="-+")
p.add_argument('-i', "--file", nargs='+',
               help='list to be converted files here')

args = p.parse_args()

for filename in args.file:
    with open(filename, "r") as csvfile:
        # Skip header:
        csvfile.readline()
        # Load to pandas:
        df = pd.read_csv(csvfile, sep="\t", header=None)
        df[3] = df[0]
        df = df.drop(columns=[4,5,8])
        df.to_csv(filename[:-4]+".bedpe", header=None, index=False, sep="\t")


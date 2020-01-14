import pandas as pd
import cooler
# import argparse


# p = argparse.ArgumentParser(description='Convert to bedpe',
#                             prefix_chars="-+")
# p.add_argument('-i', "--file", nargs='+',
#                help='list to be converted files here')

# args = p.parse_args()

# for filename in args.file:

filename = "OnTAD_cis_chr1.tad"
filep = "testdata/G2.fc_1_2.wOldG2.cis.1000.mcool"
binsize = 50000
cooler_obj = cooler.Cooler(f'{filep}::/resolutions/{binsize}')
bins = cooler_obj.bins()[:]
pixels = cooler_obj.pixels()[:]
csvfile = open(filename, "r")
df = pd.read_csv(csvfile, sep="\t", header=None)
type(pixels)
print(pixels)
print(df)
df = df.rename(columns={0: "bin1_id", 1: "bin2_id" })
print(df)
type(df)
result = cooler.annotate(df, bins)
print(result)
df=df.drop(columns=["weight1","weight2"])


with open(filename, "r") as csvfile:
    # Skip header:
    # csvfile.readline()
    # Load to pandas:
    df = pd.read_csv(csvfile, sep="\t", header=None)
    # duplicate coordinates
    df[5] = df[0]
    df[6] = df[1]
    # reorder for badpe format
    df = df[[0, 1, 5, 6, 2, 3, 4]]
    df.to_csv(filename[:-4]+".bedpe", header=None, index=False, sep="\t")



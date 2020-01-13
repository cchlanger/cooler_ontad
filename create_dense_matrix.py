"""Command line interface for OnTAD preprocessing"""
import click
import os
import cooler
# import cooltools
# import cooltools.snipping
import bioframe
from NGS import HiCTools as HT
import numpy as np


# set wd
os.chdir(".")

# get chromosomal arms

chromsizes = bioframe.fetch_chromsizes("hg19")
arms = HT.getArmsHg19()

# load in data

# fileP = "testdata/G2.fc_1_2.wOldG2.cis.1000.mcool"
# binsize = 50000


@click.command()
@click.argument('filep', type=click.Path(exists=True))
@click.option('--binsize', '-b', 'binsize', default=50000,
              help='Resulution size.',
              show_default=True)
def create_dense_matrix(filep, binsize):
    cooler_obj = cooler.Cooler(f'{filep}::/resolutions/{binsize}')
    filename = os.path.basename(filep)
    filename_base = filename[:-6]
    for chr_name in chromsizes.index:
        matrix = cooler_obj.matrix(balance=True).fetch(chr_name)
        matrix = np.nan_to_num(matrix)
        print(chr_name)
        print("tmpdata/%s.%s.matrix" % (filename_base, chr_name))
        np.savetxt("tmpdata/%s.%s.matrix" % (filename_base, chr_name), matrix, fmt='%1.8f', delimiter="\t")

if __name__ == "__main__":
    create_dense_matrix()

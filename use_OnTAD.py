"""Command line interface for OnTAD preprocessing"""
import click
import os
import cooler
import tempfile
# import cooltools
# import cooltools.snipping
import bioframe
from NGS import HiCTools as HT
import numpy as np
import pandas as pd

import subprocess
# import logging
# import shlex
import glob
import re
import shutil

# get chromosomal arms

chromsizes = bioframe.fetch_chromsizes("hg19")
arms = HT.getArmsHg19()


def convert_to_bedpe(cooler_filename, binsize, tad_folder):
    # get bins from cooler file
    binsize = 50000
    cooler_obj = cooler.Cooler(f'{cooler_filename}::/resolutions/{binsize}')
    bins = cooler_obj.bins()[:]
    f_list = glob.glob('%s/*.tad' % tad_folder)
    # pixels = cooler_obj.pixels()[:]
    results = []
    for tad_filename in f_list:
        if os.stat(tad_filename).st_size != 0:
            with open(tad_filename, "r") as csvfile:
                csvfile = open(tad_filename, "r")
                df = pd.read_csv(csvfile, sep="\t", header=None)
                df = df.rename(columns={0: "bin1_id", 1: "bin2_id"})
                # Substract one from all bin1_id and bin2_id to correct for one-based indexing in OnTAD
                df = df.sub([1, 1, 0, 0, 0], axis="columns")
                # Extract Chromosome name out of filename
                region = re.findall(r"(chr\d|chr\w)", tad_filename)[0]
                # Add offset for particular chromosome
                offset = cooler_obj.offset(region)
                df = df.add([offset, offset, 0, 0, 0], axis="columns")
                # Use cooler anotate
                result = cooler.annotate(df, bins)
                result = result.drop(columns=["weight1", "weight2"])
                results.append(result)
                # Output bedpde with addtitonal columns TADlevel  TADmean  TADscore
    bedpe = pd.concat(results)
    bedpe.to_csv(os.path.basename(cooler_filename)[:-6] + ".binsize_" + str(binsize) + ".bedpe", header=None, index=False, sep="\t")

# TODO:Stitch files together


def create_dense_matrix(filep, binsize):
    cooler_obj = cooler.Cooler(f'{filep}::/resolutions/{binsize}')
    filename = os.path.basename(filep)
    filename_base = filename[:-6]
    temp_folder = tempfile.mkdtemp(suffix=None, prefix=None, dir=None)
    for chr_name in chromsizes.index:
        matrix = cooler_obj.matrix(balance=True).fetch(chr_name)
        matrix = np.nan_to_num(matrix)
        # print(chr_name)
        # print("tmpdata/%s.%s.matrix" % (temp_folder ,filename_base, chr_name))
        np.savetxt("%s/%s.%s.matrix" % (temp_folder, filename_base, chr_name), matrix, fmt='%1.8f', delimiter="\t")
    return temp_folder


@click.command()
@click.argument('filep', type=click.Path(exists=True))
@click.option('--binsize', default=50000,
              help='Resulution size.',
              show_default=True)
@click.option('--penalty', default=0.1,
              help='-penalty <float> The penalty applied in scoring function to select positive TADs. Higher penalty score will result in fewer TADs.',
              show_default=True)
@click.option('--maxsz', default=200,
              help='-maxsz <int> The maximum size of TADs can be called. The size is determined by number of bins covered in the contact matrix.',
              show_default=True)
@click.option('--minsz', default=3,
              help='-minsz <int> The minimum size of TADs can be called. The size is determined by number of bins covered in the contact matrix.',
              show_default=True)
@click.option('--ldiff', default=1.96,
              help='-ldiff <float> The cut-off to determine local minimum. (local maximum - local minimum >= ldiff*std)',
              show_default=True)
@click.option('--lsize', default=5,
              help='-lsize <int> The local region size that used to determine local minimum',
              show_default=True)
# -log2 <boolean> if specified, log2(contact frequency) will be used to call TADs.
# -o <file path> The file path for the TAD calling results.
def main(filep, binsize, penalty, minsz, maxsz, ldiff, lsize):
    # Create Matrix for every chromosome for OnTAD
    matrix_folder = create_dense_matrix(filep, binsize)
    # Call On TAD
    f_list = glob.glob('%s/*.matrix' % matrix_folder)
    processes = []
        # Creates the files in the format for OnTADO
    tad_folder = tempfile.mkdtemp(suffix=None, prefix=None, dir=None)
    for filename in f_list:
        # TODO parse all options
        filename_base = os.path.basename(filename)[:-7]
        command_string = "OnTAD %s -penalty 0.1 -maxsz 200 -o %s/%s" % (filename, tad_folder, filename_base)
        processes.append(subprocess.Popen(command_string, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True))
    for proc in processes:
        output, error = proc.communicate()
        if proc.returncode != 0:
            print(error.decode("utf-8"))
    # Creates the .bedpe wfor all chromosomes
    convert_to_bedpe(filep, binsize, tad_folder)
    # remove /temp dirs
    shutil.rmtree(matrix_folder)
    shutil.rmtree(tad_folder)


if __name__ == "__main__":
    main()

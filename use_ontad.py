"""Command line interface for OnTAD preprocessing"""
import click
import os
import cooler
import tempfile
import bioframe
import numpy as np
import pandas as pd
import subprocess
import glob
import re
import shutil
import warnings

# get chromosomal arms

chromsizes = bioframe.fetch_chromsizes("hg19")


def convert_to_bedpe(cooler_filename, binsize, tad_folder, penalty, minsz, maxsz, ldiff, lsize, output, short_name):
    # get bins from cooler file
    cooler_obj = cooler.Cooler(f'{cooler_filename}::/resolutions/{binsize}')
    bins = cooler_obj.bins()[:]
    f_list = glob.glob('%s/*.tad' % tad_folder)
    # pixels = cooler_obj.pixels()[:]
    results = []
    for tad_filename in f_list:
        if os.stat(tad_filename).st_size != 0:
            with open(tad_filename, "r") as csvfile:
                csvfile = open(tad_filename, "r")
                # Skip first row, since this encodes level 0, which is the whole chromosome and check whether the file is empty afterwards
                try:
                    df = pd.read_csv(csvfile, sep="\t", header=None, skiprows=1)
                except pd.errors.EmptyDataError:
                    warnings.warn(f"{tad_filename} contains no TADs!")
                    continue
                df = df.rename(columns={0: "bin1_id", 1: "bin2_id"})
                # Substract one from all bin1_id and bin2_id to correct for one-based indexing in OnTAD
                df = df.sub([1, 1, 0, 0, 0], axis="columns")
                # Extract Chromosome name out of filename
                region = re.findall(r"(chr\d+|chr\w)", tad_filename)[0]
                # Add offset for particular chromosome
                offset = cooler_obj.offset(region)
                df = df.add([offset, offset, 0, 0, 0], axis="columns")
                # Use cooler anotate
                result = cooler.annotate(df, bins)
                result = result.drop(columns=["weight1", "weight2", "bin1_id", "bin2_id"])
                # For bedpe we need mitpoint of the bins and both points are identical
                mid_pos1 = (result["start1"] + result["end1"]) / 2
                mid_pos2 = (result["start2"] + result["end2"]) / 2
                result["start1"] = mid_pos1.astype('int32')
                result["end1"] = mid_pos2.astype('int32')
                result["start2"] = mid_pos1.astype('int32')
                result["end2"] = mid_pos2.astype('int32')
                results.append(result)
                # Output bedpde with addtitonal columns TADlevel  TADmean  TADscore
    bedpe = pd.concat(results)
    # Only adds paremeters to basefilename if not default or specified by the user
    string_binsize = ""
    string_penalty = ""
    string_minsz = ""
    string_maxsz = ""
    string_ldiff = ""
    string_lsize = ""
    string_output = ""
    if (binsize != 50000 or short_name is not True):
        string_binsize = ".binsize_" + str(binsize)
    if (penalty != 0.1 or short_name is not True):
        string_penalty = ".penalty_" + str(penalty)
    if (minsz != 3 or short_name is not True):
        string_minsz = ".minsz_" + str(minsz)
    if (maxsz != 200 or short_name is not True):
        string_maxsz = ".maxsz_" + str(maxsz)
    if (ldiff != 1.96 or short_name is not True):
        string_ldiff = ".ldiff_" + str(ldiff)
    if (lsize != 5 or short_name is not True):
        string_lsize = ".lsize_" + str(lsize)
    if output is not None:
        bedpe.to_csv(output, header=None, index=False, sep="\t")
    else:
        bedpe.to_csv(os.path.basename(cooler_filename)[:-6] + string_binsize + string_penalty + string_minsz + string_maxsz + string_ldiff + string_lsize + string_output + ".bedpe", header=None, index=False, sep="\t")


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
              help='--penalty <float> The penalty applied in scoring function to select positive TADs. Higher penalty score will result in fewer TADs.',
              show_default=True)
@click.option('--maxsz', default=200,
              help='--maxsz <int> The maximum size of TADs can be called. The size is determined by number of bins covered in the contact matrix.',
              show_default=True)
@click.option('--minsz', default=3,
              help='--minsz <int> The minimum size of TADs can be called. The size is determined by number of bins covered in the contact matrix.',
              show_default=True)
@click.option('--ldiff', default=1.96,
              help='--ldiff <float> The cut-off to determine local minimum. (local maximum - local minimum >= ldiff*std)',
              show_default=True)
@click.option('--lsize', default=5,
              help='--lsize <int> The local region size that used to determine local minimum',
              show_default=True)
@click.option('--dense_matrix_only', is_flag=True,
              help='--dense_matrix_only If chosen a folder dense_matrices is created and further processing is stoped')
@click.option('--o', default=None,
              help='--o <strig> Output filename. If left empty .mcool basename is extened with all non default parameters.')
@click.option('--short_name', is_flag=True,
              help='--short_name If chosen only paremeters are added to basefilename if they are not default or specified by the user.')
# -log2 <boolean> if specified, log2(contact frequency) will be used to call TADs.
# -o <file path> The file path for the TAD calling results.
def main(filep, binsize, penalty, minsz, maxsz, ldiff, lsize, dense_matrix_only, o, short_name):
    # Create Matrix for every chromosome for OnTAD
    print("Creating Matrix ...")
    matrix_folder = create_dense_matrix(filep, binsize)
    # Call OnTAD
    if (dense_matrix_only is True):
        print("Saving Matrix then Aborting!")
        if o is None:
            shutil.copytree(matrix_folder, "dense_matrices")
        else:
            shutil.copytree(matrix_folder, os.path.dirname(o) + "/dense_matrices")
        return 0
    f_list = glob.glob('%s/*.matrix' % matrix_folder)
    processes = []
    # Runs OnTAD
    print("Running OnTAD ...")
    tad_folder = tempfile.mkdtemp(suffix=None, prefix=None, dir=None)
    for filename in f_list:
        filename_base = os.path.basename(filename)[:-7]
        command_string = "OnTAD %s -penalty %s -maxsz %s -maxsz %s -ldiff %s -lsize %s -o %s/%s" % (filename, penalty, minsz, maxsz, ldiff, lsize, tad_folder, filename_base, short_name)
        processes.append(subprocess.Popen(command_string, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True))
    for proc in processes:
        output, error = proc.communicate()
        if proc.returncode != 0:
            print(error.decode("utf-8"))
    # Creates the .bedpe for all chromosomes
    print("Creating OnTAD ...")
    convert_to_bedpe(filep, binsize, tad_folder, penalty, minsz, maxsz, ldiff, lsize, o, short_name)
    # remove /temp dirs
    # print(matrix_folder)
    # print(tad_folder)
    shutil.rmtree(matrix_folder)
    shutil.rmtree(tad_folder)


if __name__ == "__main__":
    main()

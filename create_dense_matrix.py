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

chromsizes = bioframe.fetch_chromsizes('hg19')
arms = HT.getArmsHg19()

# load in data

HICPATH = "/groups/gerlich/experiments/Experiments_004800/004812/Sequencing_data/Pooled_FC_1_2/cooler/"

BINSIZE = 50000
BARCODE = "G2.fc_1_2.wOldG2"
clrs = {interType:
        cooler.Cooler(
            os.path.join(HICPATH, f'{BARCODE}.{interType}.1000.mcool::/resolutions/{BINSIZE}'))
        for interType in ["cis", "trans"]}


cis_cooler = clrs['cis']
cis_matrix = cis_cooler.matrix(balance=True).fetch("chr18")
cis_matrix = np.nan_to_num(cis_matrix)
np.savetxt("tmpdata/G2.fc_1_2.wOldG2.chr18.cis.matrix", cis_matrix, fmt='%1.4f', delimiter="\t")

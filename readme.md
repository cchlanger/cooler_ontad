OnTAD repo:
https://github.com/anlin00007/OnTAD

Create Input with:
python create_dense_matrix.py testdata/G2.fc_1_2.wOldG2.cis.1000.mcool -b 50000

Use our docker:
singularity run docker://gerlichlab/mmhic:release-1.0

Run there tool with following parameters:
OnTAD tmpdata/G2.fc_1_2.wOldG2.cis.1000.chr1.matrix -penalty 0.1 -maxsz 200 -o OnTAD_cis_chr1



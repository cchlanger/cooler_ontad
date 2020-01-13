 OnTAD repo:
 https://github.com/anlin00007/OnTAD

Create Input with:
python create_dense_matrix.py testdata/G2.fc_1_2.wOldG2.cis.1000.mcool -b 50000


 Use there docker:
 singularity run docker://anlin00007/ontad:v1.2

 Run there tool with following parameters:
 
 OnTAD tmpdata/G2.fc_1_2.wOldG2.chr18.cis.matrix -penalty 0.1 -maxsz 200 -bedout 18 78077248 50000 -o OnTAD_cis_chr18

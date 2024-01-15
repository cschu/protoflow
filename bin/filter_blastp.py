#!/usr/bin/env python

import argparse
import csv

import pandas as pd

"""
NODE_454_length_15407_cov_11.9692_ID_907_778_2124_-     k141_929_3      100.000 448     0       0       1       448     1       448     0.0     917     448     449     100     448
NODE_454_length_15407_cov_11.9692_ID_907_778_2124_-     k141_1341_4     24.528  106     58      5       240     338     195     285     6.2     25.4    448     300     22      45
NODE_454_length_15407_cov_11.9692_ID_907_778_2124_-     k141_1388_6     20.988  81      63      1       120     199     81      161     9.2     24.6    448     172     18      33
NODE_454_length_15407_cov_11.9692_ID_907_2236_3585_-    k141_929_4      100.000 449     0       0       1       449     1       449     0.0     931     449     450     100     449
NODE_454_length_15407_cov_11.9692_ID_907_2236_3585_-    k141_1388_4     25.455  110     71      4       141     243     44      149     1.7     27.3    449     329     23      47
NODE_454_length_15407_cov_11.9692_ID_907_2236_3585_-    k141_1093_39    28.358  67      35      2       256     309     476     542     1.9     27.3    449     584     12      31
NODE_454_length_15407_cov_11.9692_ID_907_2236_3585_-    k141_1297_12    52.381  21      10      0       320     340     40      60      2.2     26.9    449     326     5       14
NODE_454_length_15407_cov_11.9692_ID_907_2236_3585_-    k141_1411_59    25.532  141     91      5       203     339     307     437     2.3     26.9    449     561     31      60
NODE_454_length_15407_cov_11.9692_ID_907_2236_3585_-    k141_934_4      30.000  60      38      2       4       63      158     213     3.0     26.6    449     318     13      31
NODE_454_length_15407_cov_11.9692_ID_907_2236_3585_-    k141_1293_1     38.710  31      15      1       423     449     1       31      3.2     25.4    449     128     6       19
NODE_454_length_15407_cov_11.9692_ID_907_2236_3585_-    k141_614_34     33.333  66      30      3       12      71      212     269     3.6     26.2    449     375     13      34
NODE_454_length_15407_cov_11.9692_ID_907_2236_3585_-    k141_1322_22    45.455  33      17      1       417     449     153     184     4.1     26.2    449     436     7       20
NODE_454_length_15407_cov_11.9692_ID_907_2236_3585_-    k141_1028_8     28.261  46      31      1       286     329     203     248     4.2     25.8    449     298     10      24
"""

def main():

	HEADER = 'qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs positive'.split(" ")

	ap = argparse.ArgumentParser()
	ap.add_argument("blast_output", type=str)
	ap.add_argument("--pident_cutoff_hi", type=float, default=97.0)
	ap.add_argument("--pident_cutoff_lo", type=float, default=33.0)
	ap.add_argument("--qcovs_cutoff_hi", type=float, default=97)
	ap.add_argument("--out_prefix", "-o", type=str, default="blastp_filtered")
	args = ap.parse_args()

	df = pd.read_csv(args.blast_output, sep="\t", columns=HEADER)

	df_hi = df[(df["pident"] > args.pident_cutoff_hi) & df["qcovs"] > args.qcovs_cutoff_hi]
	df_hi.to_csv(f"{args.out_prefix}.hi_conf.tsv", sep="\t", index=False)

	df_lo = df[(df["pident"] > args.pident_cutoff_lo) & (df["pident"] <= args.pident_cutoff_hi) & df["qcovs"] > args.qcovs_cutoff_hi]
	df_lo.to_csv(f"{args.out_prefix}.lo_conf.tsv", sep="\t", index=False)




	...

if __name__ == "__main__":
	main()
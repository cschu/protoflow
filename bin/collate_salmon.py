#!/usr/bin/env python

import argparse
import gzip

import pandas as pd


def parse_attrib(s):
	return dict(item.split("=") for item in s.strip().strip(";").split(";"))

def read_proteome(f):
	def parse_protein(line):
		# >k141_0_1 # 2 # 313 # 1 # ID=1_1;partial=11;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.577
		line = line.strip()[1:].split(" ")
		d = parse_attrib(line[8])
		# return line[8].split(";")[0].split("=")[1], line[0]
		return d.get("ID"), (line[0], d.get("partial"))

	open_f = gzip.open if f.endswith(".gz") else open
	with open_f(f, "rt") as _in:
		# return [line.strip()[1:].split(" ")[0] for line in _in if line[0] == ">"]
		return dict(parse_protein(line) for line in _in if line[0] == ">")

def read_metaGT_profiles(protein_coding_genes, metaG_profiles=None, metaT_profiles=None):
	profiles_df = pd.DataFrame(index=sorted(protein_coding_genes))
	profiles_df["metaG"] = 0.0
	profiles_df["metaG_tpm"] = 0.0
	profiles_df["metaT"] = 0.0
	profiles_df["metaT_tpm"] = 0.0
	# profiles_df["Length"] = None

	usecols = ["Name", "EffectiveLength", "NumReads"]
	has_length = False


	if metaG_profiles:
		for f in metaG_profiles:
			df = profiles_df.merge(
				pd.read_csv(f, sep="\t", index_col=0, usecols=usecols),
				left_index=True,
				right_index=True,
				how="outer",
			)
			print(df.head())
			if not has_length:
				profiles_df["Length"] = df["EffectiveLength"]
				# usecols = ["Name", "NumReads"]
			profiles_df["metaG"] += df["NumReads"]
			profiles_df["metaG_tpm"] += df["NumReads"] / (df["EffectiveLength"] / 1000)

	if metaT_profiles:
		for f in metaT_profiles:
			df = profiles_df.merge(
				pd.read_csv(f, sep="\t", index_col=0, usecols=usecols),
				left_index=True,
				right_index=True,
				how="outer",
			)
			if not has_length:
				profiles_df["Length"] = df["EffectiveLength"]
				# usecols = ["Name", "NumReads"]
			profiles_df["metaT"] += df["NumReads"]
			profiles_df["metaT_tpm"] += df["NumReads"] / (df["EffectiveLength"] / 1000)

	# https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
	# Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
	# Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
	# Divide the RPK values by the “per million” scaling factor. This gives you TPM.


	metaG_scaling = profiles_df["metaG_tpm"].sum(skipna=True, numeric_only=True) / 1e6
	metaT_scaling = profiles_df["metaT_tpm"].sum(skipna=True, numeric_only=True) / 1e6

	profiles_df["metaG_tpm"] /= metaG_scaling
	profiles_df["metaT_tpm"] /= metaT_scaling
	
	return profiles_df


def main():
	ap = argparse.ArgumentParser()
	# ap.add_argument("sample_id", type=str)
	ap.add_argument("predicted_proteins", type=str)
	ap.add_argument("--metaG_counts", nargs="*", type=str)
	ap.add_argument("--metaT_counts", nargs="*", type=str)
	ap.add_argument("--output_prefix", "-o", type=str, default="collated")
	args = ap.parse_args()

	proteome_d = read_proteome(args.predicted_proteins)

	metaGT_profiles_df = read_metaGT_profiles([pid for pid, _ in proteome_d.values()], metaG_profiles=args.metaG_counts, metaT_profiles=args.metaT_counts)

	metaGT_profiles_df.to_csv(f"{args.output_prefix}.salmon.tsv", sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
	main()
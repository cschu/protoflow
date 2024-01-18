#!/usr/bin/env python

import argparse
import gzip

import pandas as pd

def read_proteome(f):
	def parse_protein(line):
		line = line.strip()[1:].split(" ")
		return line[8].split(";")[0].split("=")[1], line[0]

	open_f = gzip.open if f.endswith(".gz") else open
	with open_f(f, "rt") as _in:
		# return [line.strip()[1:].split(" ")[0] for line in _in if line[0] == ">"]
		return dict(parse_protein(line) for line in _in if line[0] == ">")
	
def read_metaP_data(f):
	open_f = gzip.open if f.endswith(".gz") else open
	with open_f(f, "rt") as _in:
		metaP = dict([(line.strip()[1:], len(next(_in).strip())) for line in _in])
		return metaP, pd.DataFrame(data=metaP.keys(), columns=["qaccver"]).set_index("qaccver")
	
def read_blastp(f, metaP_df, hi_conf_pident_cutoff=97, lo_conf_pident_cutoff=33, qcovs_cutoff=97):
	header = 'qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs positive ppos'.split(" ")
	blastp_df = pd.read_csv(f, sep="\t", names=header)

	df_hi = blastp_df[(blastp_df["pident"] > hi_conf_pident_cutoff) & (blastp_df["qcovs"] > qcovs_cutoff)]
	df_hi = pd.merge(df_hi, metaP_df, right_on="qaccver", left_on="qaccver", how="outer").drop(["mismatch", "gapopen", "qstart", "qend", "sstart", "send"], axis=1)

	unseen_df = pd.DataFrame(df_hi[df_hi["saccver"].isna()]["qaccver"]).set_index(["qaccver"])

	df_lo = blastp_df[((blastp_df["pident"] > lo_conf_pident_cutoff) & (blastp_df["pident"] <= hi_conf_pident_cutoff)) & (blastp_df["qcovs"] > qcovs_cutoff)]
	df_lo = pd.merge(df_lo, unseen_df, right_on="qaccver", left_on="qaccver", how="inner").drop(["mismatch", "gapopen", "qstart", "qend", "sstart", "send"], axis=1)

	df_hi_clean = df_hi[df_hi["saccver"].notna()].drop_duplicates(["qaccver"], keep="first")
	df_hi_clean["confidence"] = "high"
	df_lo_clean = df_lo[df_lo["saccver"].notna()].drop_duplicates(["qaccver"], keep="first")
	df_lo_clean["confidence"] = "low"

	blastp_combined_df = pd.concat([df_hi_clean, df_lo_clean]).sort_values(["qaccver", "confidence"])

	return blastp_combined_df

def read_miniprot(f, metaP, id_proteome):

	def parse_attrib(s):
		return dict(item.split("=") for item in s.strip(";").split(";"))
	
	with open(f, "rt") as _in:
		records = []
		for line in _in:
			line = line.strip().split("\t")
			m_contig, _, m_type, m_start, m_end, _, m_strand, _, m_attrib = line[:9]
			if m_type != "mRNA":
				continue
			_, _, _, p_start, p_end, _, p_strand, _, p_attrib, overlap = line[9:]
			
			if overlap == "0":
				continue
			
			m_attrib, p_attrib = parse_attrib(m_attrib), parse_attrib(p_attrib)
			
			overlap = int(overlap)
			plen = int(metaP.get(m_attrib.get("Target")))
			pid = id_proteome.get(p_attrib.get("ID"))
			mlen = int(m_end) - int(m_start) + 1
			pplen = int(p_end) - int(p_start) + 1
			
			records.append({
				"qaccver": m_attrib.get("Target"),
				"saccver": pid,
				"pident": m_attrib.get("Identity"),
				"length": plen,
				"qcov": overlap/mlen,
				"scov": overlap/pplen,
				"positive": m_attrib.get("Identity"),
				"partial": p_attrib.get("partial", "00") != "00",
				"qlen": plen,
				"slen": pplen / 3,
				"rank": int(m_attrib.get("Rank", -1)),
			})
	
	return pd.DataFrame.from_records(records)

def read_metaGT_profiles(protein_coding_genes, metaG_profiles=None, metaT_profiles=None):
	profiles_df = pd.DataFrame(index=sorted(protein_coding_genes))
	profiles_df["metaG"] = 0.0
	profiles_df["metaT"] = 0.0
	profiles_df["Length"] = None

	if metaG_profiles:
		for f in metaG_profiles:
			df = profiles_df.merge(
				pd.read_csv(f, sep="\t", index_col=0, usecols=["Name", "Length", "NumReads"]),
				left_index=True,
				right_index=True,
				how="outer",
			)
			profiles_df["Length"] = df["Length"]
			profiles_df["metaG"] += df["NumReads"]

	if metaT_profiles:
		for f in metaT_profiles:
			df = profiles_df.merge(
				pd.read_csv(f, sep="\t", index_col=0, usecols=["Name", "Length", "NumReads"]),
				left_index=True,
				right_index=True,
				how="outer",
			)
			profiles_df["Length"] = df["Length"]
			profiles_df["metaT"] += df["NumReads"]
	
	return profiles_df
	

			










def main():
	ap = argparse.ArgumentParser()
	ap.add_argument("sample_id", type=str)
	ap.add_argument("predicted_proteins", type=str)
	ap.add_argument("metaP_proteins", type=str)
	ap.add_argument("blastp_output", type=str)
	ap.add_argument("miniprot_output", type=str)
	ap.add_argument("--metaG_counts", nargs="*", type=str)
	ap.add_argument("--metaT_counts", nargs="*", type=str)
	ap.add_argument("--output_prefix", "-o", type=str, default="collated")
	args = ap.parse_args()

	proteome_d = read_proteome(args.predicted_proteins)
	metaP_d, metaP_df = read_metaP_data(args.metaP_proteins)
	blastp_df = read_blastp(args.blast_output, metaP_df)
	miniprot_df = read_miniprot(args.miniprot_output, metaP_d, proteome_d)

	metaGT_profiles_df = read_metaGT_profiles(proteome_d.values(), metaG_profiles=args.metaG_counts, metaT_profiles=args.metaT_counts)

	metaP_df.to_csv(f"{args.output_prefix}.metaP_df.tsv", sep="\t", index=False)
	blastp_df.to_csv(f"{args.output_prefix}.blastp_df.tsv", sep="\t", index=False)
	miniprot_df.to_csv(f"{args.output_prefix}.miniprot_df.tsv", sep="\t", index=False)
	metaGT_profiles_df.to_csv(f"{args.output_prefix}.metaGT_profiles_df.tsv", sep="\t", index=False)


if __name__ == "__main__":
	main()
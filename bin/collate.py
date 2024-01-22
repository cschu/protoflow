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

def filter_miniprot(df, metaP_df, hi_conf_pident_cutoff=97, lo_conf_pident_cutoff=33, qcovs_cutoff=97):
	df_hi = df[(df["pident"] * 100 > hi_conf_pident_cutoff) & (df["qcov"] * 100 > qcovs_cutoff)]
	df_hi = pd.merge(df_hi, metaP_df, right_on="qaccver", left_on="qaccver", how="outer")  # .drop(["mismatch", "gapopen", "qstart", "qend", "sstart", "send"], axis=1)
	df_hi["confidence"] = "high"

	unseen_df = pd.DataFrame(df_hi[df_hi["saccver"].isna()]["qaccver"]).set_index(["qaccver"])
	df_hi = df_hi[df_hi["saccver"].notna()]

	df_lo = df[((df["pident"] * 100 > lo_conf_pident_cutoff) & (df["pident"] * 100 <= hi_conf_pident_cutoff)) & (df["qcov"] * 100 > qcovs_cutoff)]
	df_lo = pd.merge(df_lo, unseen_df, right_on="qaccver", left_on="qaccver", how="inner")  # .drop(["mismatch", "gapopen", "qstart", "qend", "sstart", "send"], axis=1)
	df_lo["confidence"] = "low"

	return pd.concat([df_hi, df_lo]).sort_values(["qaccver", "confidence"])

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
			

			# if float(m_attrib.get("Identity", 0)) * 100 < 97:
			# 	continue

			records.append({
				"qaccver": m_attrib.get("Target"),
				"saccver": pid,
				"pident": float(m_attrib.get("Identity")),
				"length": plen,
				"qcov": overlap/mlen,
				"scov": overlap/pplen,
				"positive": float(m_attrib.get("Positive")),
				"prodigal_partial": p_attrib.get("partial", "00") != "00",
				"qlen": plen,
				"slen": pplen / 3,
				"rank": int(m_attrib.get("Rank", -1)),
			})
	
	return pd.DataFrame.from_records(records)





def read_metaGT_profiles(protein_coding_genes, metaG_profiles=None, metaT_profiles=None):
	profiles_df = pd.DataFrame(index=sorted(protein_coding_genes))
	profiles_df["metaG"] = 0.0
	profiles_df["metaT"] = 0.0
	# profiles_df["Length"] = None

	usecols = ["Name", "Length", "NumReads"]

	if metaG_profiles:
		for f in metaG_profiles:
			df = profiles_df.merge(
				pd.read_csv(f, sep="\t", index_col=0, usecols=usecols),
				left_index=True,
				right_index=True,
				how="outer",
			)
			print(df.head())
			if len(usecols) == 3:
				profiles_df["Length"] = df["Length"]
				usecols = ["Name", "NumReads"]
			profiles_df["metaG"] += df["NumReads"]

	if metaT_profiles:
		for f in metaT_profiles:
			df = profiles_df.merge(
				pd.read_csv(f, sep="\t", index_col=0, usecols=usecols),
				left_index=True,
				right_index=True,
				how="outer",
			)
			if len(usecols) == 3:
				profiles_df["Length"] = df["Length"]
				usecols = ["Name", "NumReads"]
			profiles_df["metaT"] += df["NumReads"]
	
	return profiles_df
	

			










def main():
	ap = argparse.ArgumentParser()
	# ap.add_argument("sample_id", type=str)
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
	blastp_df = read_blastp(args.blastp_output, metaP_df)
	miniprot_df = read_miniprot(args.miniprot_output, metaP_d, proteome_d)
	miniprot_df = filter_miniprot(miniprot_df, metaP_df)

	metaGT_profiles_df = read_metaGT_profiles(proteome_d.values(), metaG_profiles=args.metaG_counts, metaT_profiles=args.metaT_counts)

	# metaP_df.to_csv(f"{args.output_prefix}.metaP_df.tsv", sep="\t", index=False)
	blastp_df.to_csv(f"{args.output_prefix}.blastp_df.tsv", sep="\t", index=False)
	miniprot_df.to_csv(f"{args.output_prefix}.miniprot_df.tsv", sep="\t", index=False)
	metaGT_profiles_df.to_csv(f"{args.output_prefix}.metaGT_profiles_df.tsv", sep="\t", index=False)


	prot_combined_df = pd.merge(blastp_df, miniprot_df, on=["qaccver"], how="outer", suffixes=("_blastp", "_miniprot"))
	prot_combined_filter = (prot_combined_df["saccver_blastp"] == prot_combined_df["saccver_miniprot"]) | (prot_combined_df["saccver_blastp"].isna()) | (prot_combined_df["saccver_miniprot"].isna())

	# evidence_both = prot_combined_df[(prot_combined_df["saccver_x"] == prot_combined_df["saccver_y"]) | (prot_combined_df["saccver_x"].isna()) | (prot_combined_df["saccver_y"].isna())]
	evidence_both_df = prot_combined_df[prot_combined_filter]
	evidence_both_df = evidence_both_df.assign(prodigal_protein=evidence_both_df.saccver_blastp, qlen_aa=evidence_both_df.qlen_blastp)
	evidence_both_df["prodigal_protein"] = evidence_both_df["prodigal_protein"].fillna(evidence_both_df["saccver_miniprot"])
	evidence_both_df["qlen_aa"] = evidence_both_df["qlen_aa"].fillna(evidence_both_df["qlen_miniprot"])

	evidence_both_df = pd.merge(
		# prot_combined_df[prot_combined_filter],
		evidence_both_df,
		metaGT_profiles_df.drop(["Length"], axis=1),
		left_on=["prodigal_protein"],
		right_index=True,
		how="inner",
	).rename(
		columns={
			"qaccver": "metaP_protein",			
		}
	).drop(
	 	["saccver_blastp", "saccver_miniprot", "qlen_blastp", "qlen_miniprot",],
	 	axis=1,
	)

	# metaP_protein	pident_blastp	length_blastp	evalue	bitscore	qlen_blastp	slen_blastp	qcovs	positive_blastp	ppos	confidence_blastp	pident_miniprot	length_miniprot	qcov	scov	positive_miniprot	prodigal_partial	qlen_miniprot	slen_miniprot	rank	confidence_miniprot	prodigal_protein	metaG	metaT
	# metaP_protein	pident_blastp	length_blastp	evalue	bitscore	qlen_blastp	slen_blastp	qcovs	positive_blastp	ppos		pident_miniprot	length_miniprot	qcov	scov	positive_miniprot	qlen_miniprot	slen_miniprot	rank	

	evidence_both_df = evidence_both_df[
		[
			"metaP_protein",
			"prodigal_protein",
			"prodigal_partial",
			"metaG",
			"metaT",
			"confidence_blastp",
			"confidence_miniprot",

		]
	]


	evidence_both_df.to_csv(f"{args.output_prefix}.metaP_hits.tsv", sep="\t", index=False)

	with open(f"{args.output_prefix}.unknown_metaP.txt", "wt") as _out:
		unseen = set(metaP_d).difference(evidence_both_df["metaP_protein"].to_list())
		print(*sorted(unseen), sep="\n", file=_out)

	

	
	




if __name__ == "__main__":
	main()
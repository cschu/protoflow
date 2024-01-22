#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"

include { metaT_input; metaG_input } from "./protoflow/workflows/input"
include { salmon_index; salmon_quant } from "./protoflow/modules/profilers/salmon"
include { miniprot; intersect_miniprot } from "./protoflow/modules/align/miniprot"
include { makeblastdb; blastp; filter_blastp } from "./protoflow/modules/align/blast"
include { collate_results; extract_unknown_proteins } from "./protoflow/modules/collate/collate"


workflow {

	metaT_input(
		Channel.fromPath(params.metaT_input_dir + "/*", type: "dir")
	)
	metaG_input(
		Channel.fromPath(params.metaG_input_dir + "/*", type: "dir")
	)

	metaP_ch = Channel.fromPath(params.metaP_input_dir + "/**.faa")
		.map { file -> 
			def meta = [:]
			meta.id = file.getParent().getName()
			return tuple(meta.id, meta, file)
		}
	metaP_ch.dump(pretty: true, tag: "metaP_ch")

	annotation_ch = Channel.fromPath(params.annotation_input_dir + "/**.gz")
		.map { file ->
			return tuple(file.getParent().getName(), file)
		}
		.groupTuple(size: 3, sort: true)
		.map { sample_id, file -> 
			def meta = [:]
			meta.id = sample_id
			return tuple(sample_id, meta, file)
		}

	annotation_ch.dump(pretty: true, tag: "annotation_ch")

	assembly_ch = Channel.fromPath(params.assembly_input_dir + "/**.fna.gz")
		.map { file -> 
			def meta = [:]
			meta.id = file.getParent().getName()
			return tuple(meta.id, meta, file)
		}		
	
	nevermore_main(
		metaT_input.out.reads.concat(metaG_input.out.reads)
			.map { sample, files -> return tuple(sample.clone(), files) }
	)

	nevermore_main.out.fastqs.dump(pretty: true, tag: "prep_fastq_ch")

	all_samples = nevermore_main.out.fastqs
		.map { sample, files ->
			return sample.id.replaceAll(/\.meta[GT](\.singles)?$/, "")			
		}
		.unique()

	all_samples.dump(pretty: true, tag: "all_samples")

	transcriptomes_ch = annotation_ch
		.combine(all_samples, by: 0)
		.map { sample_id, sample, files -> return tuple(sample, [files[1]]) }

	transcriptomes_ch.dump(pretty: true, tag: "transcriptomes_ch")

	proteomes_ch = annotation_ch
		.combine(all_samples, by: 0)
		.map { sample_id, sample, files -> return tuple(sample, [files[0]]) }
	proteomes_ch.dump(pretty: true, tag: "proteomes_ch")

	gff_ch = annotation_ch
		.combine(all_samples, by: 0)
		.map { sample_id, sample, files -> return tuple(sample_id, sample, [files[2]]) }		

	makeblastdb(proteomes_ch, "prot")

	salmon_index(transcriptomes_ch)

	salmon_quant_ch = nevermore_main.output.fastqs
		.map { sample, files -> return tuple(sample.id.replaceAll(/\.meta[GT](\.singles)?$/, ""), sample.clone(), [files].flatten()) }
		.combine(salmon_index.out.index.map { sample, files -> return tuple(sample.id, files) }, by: 0)
		.map { sample_id, sample, fastqs, index ->
			return tuple(sample.clone(), fastqs, index)
		}
	
	salmon_quant_ch.dump(pretty: true, tag: "salmon_quant_ch")

	salmon_quant(salmon_quant_ch)


	blastp_ch = metaP_ch
		.map { sample_id, sample, files -> return tuple(sample_id, [files])}
		.join(
			makeblastdb.out.db.map { sample, db -> return tuple(sample.id, db) },
			by: 0
		)
		.map { sample_id, proteins, db ->
			def meta = [:]
			meta.id = sample_id
			return tuple(meta, proteins, db)
		}

	blastp_ch.dump(pretty: true, tag: "blastp_ch")

	blastp(blastp_ch)

	// filter_blastp(blastp.out.blastp)


	genomes_ch = assembly_ch
		.combine(all_samples, by: 0)
		.map { sample_id, sample, files -> return tuple(sample_id, [files])}
	genomes_ch.dump(pretty: true, tag: "genomes_ch")
	
	miniprot_ch = metaP_ch
		.map { sample_id, sample, files -> return tuple(sample_id, [files])}
		.join(genomes_ch, by: 0)
		.map { sample_id, proteins, genome ->
			def meta = [:]
			meta.id = sample_id
			return tuple(meta, proteins, genome)
		}
	miniprot_ch.dump(pretty: true, tag: "miniprot_ch")

	miniprot(miniprot_ch)

	intersect_miniprot(
		miniprot.out.gff
			.map { sample, gff -> return tuple(sample.id, sample.clone(), gff) }
			.join(gff_ch, by: 0)
			.map { sample_id, sample_x, miniprot_gff, sample_y, annotation_gff -> return tuple(sample_x.clone(), miniprot_gff, annotation_gff) }
	)

	salmon_results_ch = salmon_quant.out.salmon_counts
		.map { sample, counts ->
			sample_id = sample.id.replaceAll(/\.meta[GT](\.singles)?$/, "")
			sample_lib_id = sample.id.replaceAll(/\.singles$/, "")
			return tuple(sample_id, counts)			
		}
		.groupTuple(by: 0, size: 4, remainder: true, sort: true)
		.map { sample_id, counts ->
			def meta = [:]
			meta.id = sample_id
			return tuple(meta, counts)
		}
	salmon_results_ch.dump(pretty: true, tag: "salmon_results_ch")

	results_ch = metaP_ch
		.map { sample_id, sample, files -> return tuple(sample, [files])}
		.join(proteomes_ch)
		.join(intersect_miniprot.out.mp_intersect)
		.join(blastp.out.blastp)
		.join(salmon_results_ch)
	results_ch.dump(pretty: true, tag: "results_ch")

	collate_results(results_ch)

	extract_unknown_proteins(
		collate_results.out.unknown_metaP
			.join(metaP_ch.map { sample_id, sample, files -> return tuple(sample, [files]) }, by: 0)
	)

	// 1235  singularity exec -B /scratch -B /g/ bedtools_latest.sif bedtools intersect -a /g/scb2/bork/data/MAGs/annotations/internal_MICROB-PREDICT/psa_megahit/prodigal/MPHU23965372ST.psa_megahit.prodigal.gff.gz -b work/86/19e7407ece45ab89080ca4c9df73ea/miniprot/17_I_106_R10/17_I_106_R10.gff -wao > test.overlap.txt
 	// 1237  singularity exec -B /scratch -B /g/ bedtools_latest.sif bedtools intersect -b /g/scb2/bork/data/MAGs/annotations/internal_MICROB-PREDICT/psa_megahit/prodigal/MPHU23965372ST.psa_megahit.prodigal.gff.gz -a work/86/19e7407ece45ab89080ca4c9df73ea/miniprot/17_I_106_R10/17_I_106_R10.gff -wao > test.overlap.txt
 	// 1239  singularity exec -B /scratch -B /g/ bedtools_latest.sif bedtools intersect -b /g/scb2/bork/data/MAGs/annotations/internal_MICROB-PREDICT/psa_megahit/prodigal/MPHU23965372ST.psa_megahit.prodigal.gff.gz -a work/86/19e7407ece45ab89080ca4c9df73ea/miniprot/17_I_106_R10/17_I_106_R10.gff -v > test.no_overlap.txt

}

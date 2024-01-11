#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"

include { metaT_input; metaG_input } from "./protoflow/workflows/input"
include { salmon_index } from "./protoflow/modules/profilers/salmon"

workflow {

	metaT_input(
		Channel.fromPath(params.metaT_input_dir + "/*", type: "dir")
	)
	metaG_input(
		Channel.fromPath(params.metaG_input_dir + "/*", type: "dir")
	)

	// metaT_ch = metaT_input.out.reads		
	// metaG_ch = metaG_input.out.reads
	genes_ch = Channel.fromPath(params.gene_input_dir + "/**.fna.gz")
		.map { file ->
			meta = [:]
			meta.id = file.getParent().getName()
			return tuple(meta, file)
		}

	genes_ch.dump(pretty: true, tag: "genes_ch")

	metaG_ch = metaG_input.out.reads.map { sample, files -> return tuple(sample.id.replaceAll(/\.metaG(\.singles)?$/, ""), sample, files) }
	metaG_ch.dump(pretty: true, tag: "metaG_ch")
	metaT_ch = metaT_input.out.reads.map { sample, files -> return tuple(sample.id.replaceAll(/\.metaT(\.singles)?$/, ""), sample, files) }
	metaT_ch.dump(pretty: true, tag: "metaT_ch")

	joined_ch = metaG_input.out.reads.map { sample, files -> return tuple(sample.id.replaceAll(/\.metaG(\.singles)?$/, ""), sample, files) }
		.join(metaT_input.out.reads.map { sample, files -> return tuple(sample.id.replaceAll(/\.metaT(\.singles)?$/, ""), sample, files) }, by: 0)
		.join(genes_ch.map  { sample, files -> return tuple(sample.id, sample, files) }, by: 0)

	joined_ch.dump(pretty: true, tag: "joined_ch")

	metaT_ch = joined_ch
		.map { sample_id, x, metaG, sample, metaT, z, genes -> return tuple(sample, metaT) }
	metaG_ch = joined_ch
		.map { sample_id, sample, metaG, y, metaT, z, genes -> return tuple(sample, metaG) }

	nevermore_main(metaT_ch.concat(metaG_ch))

	// Invalid method invocation `call` with arguments: 
	// [
	// 	17_A_007_S00, 
	// 	[id:17_A_007_S00.metaG.singles, is_paired:false, library:paired, library_source:metaG], 
	// 	[/scratch/schudoma/WORK/MICROB_PREDICT_proteomics.1834/76/1f98b4bb9d812f087f85394758c80e/fastq/17_A_007_S00.metaG.singles/17_A_007_S00.metaG.singles_R1.fastq.gz], 
	// 	[id:17_A_007_S00.metaT.singles, is_paired:false, library:single, library_source:metaT], 
	// 	[/scratch/schudoma/WORK/MICROB_PREDICT_proteomics.1834/3f/78822bac0e04b37287717dcd64b1b2/fastq/17_A_007_S00.metaT.singles/17_A_007_S00.metaT.singles_R1.fastq.gz], 
	// 	[id:17_A_007_S00], 
	// 	/g/scb/bork/schudoma/tasks/MICROB_PREDICT_proteomics.1834/input/prodigal/17_A_007_S00/17_A_007_S00.fna.gz
	// ]
	// (java.util.ArrayList) on _closure10 type
	salmon_index_input_ch = joined_ch.map { sample_id, x, metaG, y, metaT, sample, genes -> return tuple(sample, [genes]) }
	
	salmon_index(salmon_index_input_ch)

}

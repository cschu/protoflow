#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"

include { metaT_input; metaG_input } from "./protoflow/workflows/input"
include { salmon_index; salmon_quant } from "./protoflow/modules/profilers/salmon"

workflow {

	metaT_input(
		Channel.fromPath(params.metaT_input_dir + "/*", type: "dir")
	)
	metaG_input(
		Channel.fromPath(params.metaG_input_dir + "/*", type: "dir")
	)

	metaP_ch = Channel.fromPath(params.metaP_input_dir + "/**.faa*")
		.map { file -> 
			meta = [:]
			meta.id = file.getParent().getName()
			return tuple(meta, file)
		}

	annotation_ch = Channel.fromPath(params.gene_input_dir + "/**.f?a.gz")
		.map { file ->
			return tuple(file.getParent().getName(), file)
		}
		.groupTuple(size: 2, sort: true)
		.map { sample_id, file -> 
			meta = [:]
			meta.id = sample_id
			return tuple(sample_id, meta, file)
		}
	
	annotation_ch.dump(pretty: true, tag: "genes_ch")
	
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
	
	salmon_index(transcriptomes_ch)

	salmon_ch = nevermore_main.output.fastqs
		.map { sample, files -> 
			return tuple(sample.id.replaceAll(/(\.singles)$/, ""), files)			
		}
		.groupTuple(by: 0, size: 2, remainder: true)
		.map { sample_id, files -> return tuple(sample_id.replaceAll(/\.meta[GT]$/, ""), sample_id, [files].flatten()) }
		.combine(salmon_index.out.index.map { sample, files -> return tuple(sample.id, files) }, by: 0)
		.map { sample_id, sample_libtype_id, files -> 
			def meta = [:]
			meta.id = sample_libtype_id
			return tuple(meta, files)
		}
	
	salmon_ch.dump(pretty: true, tag: "salmon_ch")


	// Invalid method invocation `call` with arguments: 
	// [17_A_007_S00, 17_A_007_S00.metaG, [/scratch/schudoma/WORK/MICROB_PREDICT_proteomics.1834/5a/29b8455039f7691791aa0c1f00373f/no_host/17_A_007_S00.metaG/17_A_007_S00.metaG_R1.fastq.gz, /scratch/schudoma/WORK/MICROB_PREDICT_proteomics.1834/5a/29b8455039f7691791aa0c1f00373f/no_host/17_A_007_S00.metaG/17_A_007_S00.metaG_R2.fastq.gz, /scratch/schudoma/WORK/MICROB_PREDICT_proteomics.1834/6b/9e189ed60f5b8f7a101d769f8287f1/merged/17_A_007_S00.metaG.singles_R1.fastq.gz], /scratch/schudoma/WORK/MICROB_PREDICT_proteomics.1834/48/acc91e9a9c4e39abc63babe46c3438/salmon/index/17_A_007_S00] (java.util.LinkedList) on _closure12 type


	// genes_ch.dump(pretty: true, tag: "genes_ch")

	// metaG_ch = metaG_input.out.reads.map { sample, files -> return tuple(sample.id.replaceAll(/\.metaG(\.singles)?$/, ""), sample, files) }
	// metaG_ch.dump(pretty: true, tag: "metaG_ch")
	// metaT_ch = metaT_input.out.reads.map { sample, files -> return tuple(sample.id.replaceAll(/\.metaT(\.singles)?$/, ""), sample, files) }
	// metaT_ch.dump(pretty: true, tag: "metaT_ch")

	// joined_ch = metaG_input.out.reads.map { sample, files -> return tuple(sample.id.replaceAll(/\.metaG(\.singles)?$/, ""), sample, files) }
	// 	.join(metaT_input.out.reads.map { sample, files -> return tuple(sample.id.replaceAll(/\.metaT(\.singles)?$/, ""), sample, files) }, by: 0)
	// 	.join(genes_ch.map  { sample, files -> return tuple(sample.id, sample, files) }, by: 0)

	// joined_ch.dump(pretty: true, tag: "joined_ch")

	// metaT_ch = joined_ch
	// 	.map { sample_id, x, metaG, sample, metaT, z, genes -> return tuple(sample.clone(), metaT) }
	// metaG_ch = joined_ch
	// 	.map { sample_id, sample, metaG, y, metaT, z, genes -> return tuple(sample.clone(), metaG) }



}

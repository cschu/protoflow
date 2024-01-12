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

	metaP_ch = Channel.fromPath(params.metaP_input_dir + "/**.faa*")
		.map { file -> 
			meta = [:]
			meta.id = file.getParent().getName()
			return tuple(meta, file)
		}

	annotation_ch = Channel.fromPath(params.gene_input_dir + "/**.f?a.gz")
		.map { file ->
			meta = [:]
			meta.id = file.getParent().getName()
			return tuple(meta.id, file)
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

	transcriptomes_ch.dump(pretty: true, tag: "transcriptomes_ch")
		// .map { sample_id, sample, files -> }

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


	// salmon_index_input_ch = joined_ch.map { sample_id, x, metaG, y, metaT, sample, genes -> return tuple(sample, [genes]) }
	
	// salmon_index(salmon_index_input_ch)

}

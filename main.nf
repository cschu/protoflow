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
		.join(genes_ch, by: 0)

	joined_ch.dump(pretty: true, tag: "joined_ch")

	metaT_ch = joined_ch
		.map { sample_id, sample, metaG, metaT, genes -> return tuple(sample, metaT) }
	metaG_ch = joined_ch
		.map { sample_id, sample, metaG, metaT, genes -> return tuple(sample, metaG) }

	nevermore_main(metaT_ch.concat(metaG_ch))

	
	salmon_index(joined_ch.map { sample_id, x, metaG, y, metaT, genes -> return tuple(sample, genes) })

}

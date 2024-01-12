#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { nevermore_main } from "./nevermore/workflows/nevermore"

include { metaT_input; metaG_input } from "./protoflow/workflows/input"
include { salmon_index; salmon_quant } from "./protoflow/modules/profilers/salmon"
include { miniprot } from "./protoflow/modules/align/miniprot"
include { makeblastdb } from "./protoflow/modules/align/blast"

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

	annotation_ch.dump(pretty: true, tag: "genes_ch")

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

	makeblastdb(proteomes_ch, "prot")

	salmon_index(transcriptomes_ch)

	salmon_ch = nevermore_main.output.fastqs
		.map { sample, files -> return tuple(sample.id.replaceAll(/\.meta[GT](\.singles)?$/, ""), sample.clone(), [files].flatten()) }
		.combine(salmon_index.out.index.map { sample, files -> return tuple(sample.id, files) }, by: 0)
		.map { sample_id, sample, fastqs, index ->
			return tuple(sample.clone(), fastqs, index)
		}
	
	salmon_ch.dump(pretty: true, tag: "salmon_ch")

	salmon_quant(salmon_ch)


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


}

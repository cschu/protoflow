process collate_salmon {
	input:
	tuple val(sample), path(proteins), path(salmon_results)

	output:
	tuple val(sample), path("salmon/collated/${sample.id}/${sample.id}.salmon.tsv"), emit: collated

	script:
	metaG_files = salmon_results.findAll( { it.name.matches("(.*)metaG(.*)") } )
	metaT_files = salmon_results.findAll( { it.name.matches("(.*)metaT(.*)") } )

	def metaG_input = ""
	if (metaG_files.size() != 0) {
		metaG_input += "--metaG_counts "
		metaG_input += metaG_files.join(' ')
	}
	def metaT_input = ""
	if (metaT_files.size() != 0) {
		metaT_input += "--metaT_counts "
		metaT_input += metaT_files.join(' ')
	}

	"""
	mkdir -p salmon/collated/${sample.id}

	collate_salmon.py -o salmon/collated/${sample.id} ${proteins} ${metaG_input} ${metaT_input}
	"""
}

process collate_results {

	input:
	tuple val(sample), path(metaP_proteins), path(proteins), path(miniprot_results), path(blastp_results), path(salmon_results)

	output:
	tuple val(sample), path("collated/${sample.id}/${sample.id}.metaP_hits.tsv"), emit: collated
	tuple val(sample), path("collated/${sample.id}/${sample.id}.no_metaP_hits.tsv"), emit: no_metaP
	tuple val(sample), path("collated/${sample.id}/${sample.id}.unknown_metaP.txt"), emit: unknown_metaP

	script:
	metaG_files = salmon_results.findAll( { it.name.matches("(.*)metaG(.*)") } )
	metaT_files = salmon_results.findAll( { it.name.matches("(.*)metaT(.*)") } )

	def metaG_input = ""
	if (metaG_files.size() != 0) {
		metaG_input += "--metaG_counts "
		metaG_input += metaG_files.join(' ')
	}
	def metaT_input = ""
	if (metaT_files.size() != 0) {
		metaT_input += "--metaT_counts "
		metaT_input += metaT_files.join(' ')
	}    

	"""
	mkdir -p collated/${sample.id}/

	collate.py -o collated/${sample.id}/${sample.id} ${proteins} ${metaP_proteins} ${blastp_results} ${miniprot_results} ${metaG_input} ${metaT_input}
	"""
	

}

process extract_unknown_proteins {
	input:
	tuple val(sample), path(metaP_proteins), path(unknown_proteins)

	output:
	tuple val(sample), path("unknown/${sample.id}/${sample.id}.faa"), emit: unknown_proteins

	script:
	"""
	mkdir -p unknown/${sample.id}/
	seqtk subseq ${metaP_proteins} ${unknown_proteins} > unknown/${sample.id}/${sample.id}.faa
	"""

}


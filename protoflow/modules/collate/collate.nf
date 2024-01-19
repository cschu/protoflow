process collate_results {

	input:
	tuple val(sample), path(metaP_proteins), path(proteins), path(miniprot_results), path(blastp_results), path(salmon_results)

	output:
	tuple val(sample), path("collated/${sample.id}/${sample.id}.*tsv"), emit: collated

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


// process collate_protein_hits {

// 	input:
// 	tuple val(sample), path(proteomics_proteins), path(hi_conf_proteins), path(lo_conf_proteins)
	
// 	output:
// 	tuple val(sample), path("collated_proteins/${sample.id}/${sample.id}.protein_hits.tsv"), emit: protein_hits

// 	script:
// 	"""
// 	mkdir -p collated_proteins/${sample.id}/


// 	"""



// }
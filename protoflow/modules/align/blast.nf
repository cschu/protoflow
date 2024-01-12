process makeblastdb {

	input:
	tuple val(sample), path(sequences)
	val dbtype

	script:
	"""
	mkdir -p blastp/${sample.id}/

	makeblastdb -dbtype ${dbtype} -parse_seqids -in ${sequences} -out blastp/${sample.id}/${sample.id}
	"""

}
process makeblastdb {

	input:
	tuple val(sample), path(sequences)
	val dbtype

	script:
	"""
	mkdir -p blastp/${sample.id}/

	gzip -dc ${sequences} > db.faa

	makeblastdb -dbtype ${dbtype} -parse_seqids -in db.faa -out blastp/${sample.id}/${sample.id}
	"""

}
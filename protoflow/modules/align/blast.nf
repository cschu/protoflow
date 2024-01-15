process makeblastdb {

	input:
	tuple val(sample), path(sequences)
	val dbtype

	output:
	tuple val(sample), path("blast/db/${sample.id}/*"), emit: db


	script:
	"""
	mkdir -p blast/db/${sample.id}/

	gzip -dc ${sequences} > db.faa

	makeblastdb -dbtype ${dbtype} -parse_seqids -in db.faa -out blast/db/${sample.id}/${sample.id}
	"""

}

process blastp {

	input:
	tuple val(sample), path(proteins), path(db)

	output:
	tuple val(sample), path("blast/blastp/${sample.id}/${sample.id}.tsv"), emit: blastp

	script:

	def outfmt = """\"7 'qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs qcovhsp'\""""


	"""
	mkdir -p blast/blastp/${sample.id}/
	blastp -num_threads ${task.cpus} -db ${sample.id} -query ${proteins} -outfmt ${outfmt} > blast/blastp/${sample.id}/${sample.id}.tsv
	"""

}
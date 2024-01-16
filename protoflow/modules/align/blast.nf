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
	tag "${sample.id}"

	input:
	tuple val(sample), path(proteins), path(db)

	output:
	tuple val(sample), path("blast/blastp/${sample.id}/${sample.id}.tsv"), emit: blastp

	script:

	// std = 'qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore'
	def outfmt = """'6 std qlen slen qcovs positive ppos'"""


	"""
	mkdir -p blast/blastp/${sample.id}/
	blastp -num_threads ${task.cpus} -db ${sample.id} -query ${proteins} -outfmt ${outfmt} > blast/blastp/${sample.id}/${sample.id}.tsv
	"""

}

process filter_blastp {
	tag "${sample.id}"

	input:
	tuple val(sample), path(blastp_table)

	output:
	tuple val(sample), path("blast/blastp_filtered/${sample.id}/${sample.id}.hi_conf.tsv"), emit: hi_conf_proteins
	tuple val(sample), path("blast/blastp_filtered/${sample.id}/${sample.id}.lo_conf.tsv"), emit: lo_conf_proteins

	script:

	"""
	mkdir -p blast/blastp_filtered/${sample.id}/

	filter_blastp.py ${blastp_table} -o blast/blastp_filtered/${sample.id}/${sample.id}
	"""

}
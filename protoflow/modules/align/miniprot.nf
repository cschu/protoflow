process miniprot {

	input:
	tuple val(sample), path(proteins), path(genome)

	script:
	"""
	mkdir -p miniprot/${sample.id}/

	miniprot -t ${task.cpus} --gff ${genome} ${proteins} >> miniprot/${sample.id}/${sample.id}.gff
	"""

}
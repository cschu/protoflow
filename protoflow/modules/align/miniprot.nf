process miniprot {

	input:
	tuple val(sample), path(proteins), path(genome)

	output:
	tuple val(sample), path("miniprot/${sample.id}/${sample.id}.gff"), emit: gff

	script:
	"""
	mkdir -p miniprot/${sample.id}/

	miniprot -t ${task.cpus} --gff ${genome} ${proteins} >> miniprot/${sample.id}/${sample.id}.gff
	"""

}

process intersect_miniprot {
	input:
	tuple val(sample), path(protein_gff), path(genome_gff)

	output:
	tuple val(sample), path("miniprot/intersect/${sample.id}/${sample.id}.gff"), emit: mp_intersect, optional: true
	tuple val(sample), path("miniprot/intersect/${sample.id}/${sample.id}.MP_INTERSECT_DONE"), emit: processed

	script:
	"""
	set -e -o pipefail

	mkdir -p miniprot/intersect/${sample.id}/

	sed "s/ /;coords=/" ${protein_gff} | sed "s/ /-/" | grep -v "#" > proteins.gff

	if [[ -s proteins.gff ]]; then
		gzip -dc ${genome_gff} | grep -v "#" > genome.gff
		bedtools intersect -wao -a proteins.gff -b genome.gff > miniprot/intersect/${sample.id}/${sample.id}.gff
	fi

	touch miniprot/intersect/${sample.id}/${sample.id}.MP_INTERSECT_DONE
	
	"""
	// 1237  singularity exec -B /scratch -B /g/ bedtools_latest.sif bedtools intersect -b /g/scb2/bork/data/MAGs/annotations/internal_MICROB-PREDICT/psa_megahit/prodigal/MPHU23965372ST.psa_megahit.prodigal.gff.gz -a work/86/19e7407ece45ab89080ca4c9df73ea/miniprot/17_I_106_R10/17_I_106_R10.gff -wao > test.overlap.txt


}
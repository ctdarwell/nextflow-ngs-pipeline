nextflow.enable.dsl=2

params.reads = "${projectDir}/data/Hd4*_{1,2}.fastq.gz"
params.ref   = "${projectDir}/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.outdir = "results"

workflow {

    reads_ch = Channel.fromFilePairs(params.reads)
    ref_ch   = Channel.fromPath("${params.ref}*").collect()

    qc_out = FASTQC(reads_ch)

    sam_ch = ALIGN(reads_ch, ref_ch)
    bam_ch = SAMTOOLS_VIEW(sam_ch)
    stats  = FLAGSTAT(bam_ch)

    all_qc = qc_out.mix(stats)

    MULTIQC(all_qc.collect())
}

process FASTQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{html,zip}"

    script:
    """
    fastqc ${reads} -o .
    """
}

process ALIGN {
    tag "$sample_id"

    cpus 3

    input:
    tuple val(sample_id), path(reads)
    path ref_files

    output:
    path "${sample_id}.sam"

    script:
    """
    bwa mem -t ${task.cpus} ${ref_files[0]} ${reads[0]} ${reads[1]} > ${sample_id}.sam
    """
}

process SAMTOOLS_VIEW {
    publishDir "${params.outdir}/bam", mode: 'copy'
	tag "${sam.baseName}"

    input:
    path sam

    output:
    path "${sam.baseName}.bam"

    script:
    """
    samtools view -Sb ${sam} > ${sam.baseName}.bam
    """
}

process FLAGSTAT {
    publishDir "${params.outdir}/flagstat", mode: 'copy'
	tag "${bam.baseName}"

    input:
    path bam

    output:
    path "*.flagstat"

    script:
    """
    samtools flagstat ${bam} > ${bam}.flagstat
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_files

    output:
    path "multiqc_report.html"

    script:
    """
	multiqc ${qc_files} -o .
    """
}

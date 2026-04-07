nextflow.enable.dsl=2

params.reads  = "${projectDir}/data/Hd4*_{1,2}.fastq.gz"
params.ref    = "${projectDir}/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.outdir = "results"

workflow {

    // input reads
    reads_ch = Channel.fromFilePairs(params.reads)

    // bundle reference + index files
    ref_ch = Channel.fromPath("${params.ref}*").collect()

    // QC
    qc_out = FASTQC(reads_ch)

    // Alignment
    sam_ch = ALIGN(reads_ch, ref_ch)

    // BAM processing
    bam_ch     = SAMTOOLS_VIEW(sam_ch)
    sorted_ch  = SORT_BAM(bam_ch)
    indexed_ch = INDEX_BAM(sorted_ch)

    // Stats (on sorted BAM)
    stats_ch = FLAGSTAT(sorted_ch)

    // Collect QC inputs properly
    fastqc_files = qc_out.map { it[1] }
    qc_files = fastqc_files.mix(stats_ch).collect()

    // MultiQC
    MULTIQC(qc_files)
}

process FASTQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.*")

    script:
    """
    fastqc ${reads} -o .
    """
}

process ALIGN {
    tag "${sample_id}"
    cpus 3

    input:
    tuple val(sample_id), path(reads)
    path ref_files

    output:
    tuple val(sample_id), path("${sample_id}.sam")

    script:
    """
    bwa mem -t ${task.cpus} ${ref_files[0]} ${reads[0]} ${reads[1]} > ${sample_id}.sam
    """
}

process SAMTOOLS_VIEW {
    publishDir "${params.outdir}/bam", mode: 'copy'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(sam)

    output:
    tuple val(sample_id), path("${sample_id}.bam")

    script:
    """
    samtools view -Sb ${sam} > ${sample_id}.bam
    """
}

process SORT_BAM {
    publishDir "${params.outdir}/bam", mode: 'copy'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
    """
    samtools sort ${bam} -o ${sample_id}.sorted.bam
    """
}

process INDEX_BAM {
    publishDir "${params.outdir}/bam", mode: 'copy'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam.bai")

    script:
    """
    samtools index ${bam}
    """
}

process FLAGSTAT {
    publishDir "${params.outdir}/flagstat", mode: 'copy'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam)

    output:
    path "${sample_id}.flagstat"

    script:
    """
    samtools flagstat ${bam} > ${sample_id}.flagstat
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
    multiqc . -o .
    """
}

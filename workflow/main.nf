nextflow.enable.dsl=2

/*
 * Process 1: Build STAR Index
 */
process STAR_INDEX {
    tag "Indexing Reference"
    publishDir "${params.base_dir}/data/reference/star_index", mode: 'copy'

    input:
    path fasta
    path gtf

    output:
    path "star_index_dir", emit: index

    script:
    """
    mkdir star_index_dir
    STAR --runMode genomeGenerate \\
         --genomeDir star_index_dir \\
         --genomeFastaFiles ${fasta} \\
         --sjdbGTFfile ${gtf} \\
         --sjdbOverhang 99 \\
         --runThreadN ${task.cpus}
    """
}

/*
 * Main Workflow Control
 */
workflow {
    // 1. Initialize Channels
    ch_fasta = file(params.genome_fasta)
    ch_gtf   = file(params.genome_gtf)

    // 2. Parse Samplesheet: [sample_id, [fastq_1, fastq_2], condition]
    ch_samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row -> 
            def meta = [id: row.sample, condition: row.condition]
            def fastqs = [ file("${params.base_dir}/${row.fastq_1}"), file("${params.base_dir}/${row.fastq_2}") ]
            return [ meta, fastqs ]
        }

    // 3. Execute Indexing
    STAR_INDEX(ch_fasta, ch_gtf)

    // Progress Reporting
    ch_samples.view { meta, reads -> "Sample: ${meta.id} [${meta.condition}] | Reads: ${reads.join(', ')}" }
}
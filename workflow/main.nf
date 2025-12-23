nextflow.enable.dsl=2

/*
 * ============================================================================
 * RNA-seq Pipeline with Dual Quantification: STAR + Salmon
 * ============================================================================
 */

/*
 * Process 1: FastQC - Raw Quality Control
 */
process FASTQC_RAW {
    tag "${meta.id}"
    publishDir "${params.outdir}/fastqc/raw", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "*.{html,zip}", emit: fastqc_reports

    script:
    """
    fastqc -t ${task.cpus} ${reads[0]} ${reads[1]}
    """
}

/*
 * Process 2: Fastp - Trimming and Filtering
 */
process FASTP {
    tag "${meta.id}"
    publishDir "${params.outdir}/fastp", mode: 'copy', pattern: "*.{json,html}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed_R{1,2}.fastq.gz"), emit: trimmed_reads
    path "*.{json,html}", emit: fastp_reports

    script:
    """
    fastp \\
        -i ${reads[0]} \\
        -I ${reads[1]} \\
        -o ${meta.id}_trimmed_R1.fastq.gz \\
        -O ${meta.id}_trimmed_R2.fastq.gz \\
        --detect_adapter_for_pe \\
        --thread ${task.cpus} \\
        --json ${meta.id}_fastp.json \\
        --html ${meta.id}_fastp.html
    """
}

/*
 * Process 3: FastQC - Trimmed Quality Control
 */
process FASTQC_TRIMMED {
    tag "${meta.id}"
    publishDir "${params.outdir}/fastqc/trimmed", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "*.{html,zip}", emit: fastqc_reports

    script:
    """
    fastqc -t ${task.cpus} ${reads[0]} ${reads[1]}
    """
}

/*
 * ============================================================================
 * SALMON WORKFLOW (Pseudo-alignment)
 * ============================================================================
 */

/*
 * Process 4: Build Salmon Index
 */
process SALMON_INDEX {
    tag "Building Salmon Index"
    publishDir "${params.outdir}/reference/salmon_index", mode: 'copy'

    input:
    path transcriptome

    output:
    path "salmon_index", emit: index

    script:
    """
    # Extract transcripts from genome using GTF
    salmon index \\
        -t ${transcriptome} \\
        -i salmon_index \\
        -k 31 \\
        --threads ${task.cpus}
    """
}

/*
 * Process 5: Salmon Quantification
 */
process SALMON_QUANT {
    tag "${meta.id}"
    publishDir "${params.outdir}/salmon", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path salmon_index

    output:
    tuple val(meta), path("${meta.id}_salmon"), emit: quant
    path "${meta.id}_salmon/quant.sf", emit: quant_sf

    script:
    """
    salmon quant \\
        -i ${salmon_index} \\
        -l A \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --validateMappings \\
        --gcBias \\
        --numBootstraps 30 \\
        -o ${meta.id}_salmon \\
        --threads ${task.cpus}
    """
}

/*
 * ============================================================================
 * STAR WORKFLOW (Alignment-based) - OPTIONAL if index exists
 * ============================================================================
 */

/*
 * Process 6: STAR Index (Optional - skip if too memory intensive)
 */
process STAR_INDEX {
    tag "Indexing Reference"
    publishDir "${params.outdir}/reference/star_index", mode: 'copy'

    input:
    path fasta
    path gtf

    output:
    path "star_index_dir", emit: index

    when:
    params.run_star == true

    script:
    """
    mkdir star_index_dir
    STAR --runMode genomeGenerate \\
         --genomeDir star_index_dir \\
         --genomeFastaFiles ${fasta} \\
         --sjdbGTFfile ${gtf} \\
         --sjdbOverhang 99 \\
         --limitGenomeGenerateRAM 45000000000 \\
         --runThreadN ${task.cpus}
    """
}

/*
 * Process 7: STAR Alignment
 */
process STAR_ALIGN {
    tag "${meta.id}"
    publishDir "${params.outdir}/alignments", mode: 'copy', 
               pattern: "*.{bam,bam.bai,Log.final.out}"

    input:
    tuple val(meta), path(reads)
    path star_index

    output:
    tuple val(meta), path("${meta.id}_Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("${meta.id}_Aligned.sortedByCoord.out.bam.bai"), emit: bai
    path "${meta.id}_Log.final.out", emit: log

    when:
    params.run_star == true

    script:
    """
    STAR --genomeDir ${star_index} \\
         --readFilesIn ${reads[0]} ${reads[1]} \\
         --readFilesCommand zcat \\
         --outFileNamePrefix ${meta.id}_ \\
         --outSAMtype BAM SortedByCoordinate \\
         --outSAMunmapped Within \\
         --outSAMattributes Standard \\
         --limitBAMsortRAM 10000000000 \\
         --runThreadN ${task.cpus}
    
    samtools index ${meta.id}_Aligned.sortedByCoord.out.bam
    """
}

/*
 * Process 8: featureCounts
 */
process FEATURECOUNTS {
    tag "Counting all samples"
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path bams
    path gtf

    output:
    path "gene_counts.txt", emit: counts
    path "gene_counts.txt.summary", emit: summary

    when:
    params.run_star == true

    script:
    def bam_list = bams.join(' ')
    """
    featureCounts \\
        -p \\
        -T ${task.cpus} \\
        -t exon \\
        -g gene_id \\
        -a ${gtf} \\
        -o gene_counts.txt \\
        ${bam_list}
    """
}

/*
 * ============================================================================
 * MAIN WORKFLOW
 * ============================================================================
 */

workflow {
    // Initialize Reference Channels
    ch_fasta = file(params.genome_fasta)
    ch_gtf   = file(params.genome_gtf)
    ch_transcriptome = file(params.transcriptome)

    // Parse Samplesheet
    ch_samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row -> 
            def meta = [id: row.sample, condition: row.condition]
            def fastqs = [ 
                file("${params.base_dir}/${row.fastq_1}"), 
                file("${params.base_dir}/${row.fastq_2}") 
            ]
            return [ meta, fastqs ]
        }

    // ========================================================================
    // QC and Preprocessing (Both workflows)
    // ========================================================================
    FASTQC_RAW(ch_samples)
    FASTP(ch_samples)
    FASTQC_TRIMMED(FASTP.out.trimmed_reads)

    // ========================================================================
    // SALMON WORKFLOW (Always runs - fast and memory efficient)
    // ========================================================================
    SALMON_INDEX(ch_transcriptome)
    SALMON_QUANT(
        FASTP.out.trimmed_reads,
        SALMON_INDEX.out.index
    )

    // ========================================================================
    // STAR WORKFLOW (Optional - controlled by params.run_star)
    // ========================================================================
    if (params.run_star) {
        STAR_INDEX(ch_fasta, ch_gtf)
        
        STAR_ALIGN(
            FASTP.out.trimmed_reads,
            STAR_INDEX.out.index
        )

        ch_all_bams = STAR_ALIGN.out.bam
            .map { meta, bam -> bam }
            .collect()
        
        FEATURECOUNTS(
            ch_all_bams,
            ch_gtf
        )
    }

    // Progress Reporting
    ch_samples.view { meta, reads -> 
        "✓ Sample: ${meta.id} [${meta.condition}]" 
    }
}

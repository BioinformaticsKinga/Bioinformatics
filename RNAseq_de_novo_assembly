nextflow.enable.dsl=2

// Define input parameters
params.reads = '/path/to/zygmukin_lupinus_luteus_reads/*.fastq.gz'
params.outdir = './zygmukin_lupinus_luteus_results'
params.transcriptome = '/path/to/zygmukin_lupinus_luteus_trinity_output/Trinity.fasta'  // Path to Trinity assembly

// Process 1: FastQC (Quality Control)
process fastqc {
    container 'biocontainers/fastqc:v0.11.9'
    input:
    path reads

    output:
    path 'zygmukin_fastqc_reports'

    script:
    """
    fastqc -o zygmukin_fastqc_reports $reads
    """
}

// Process 2: Trimmomatic (Read Trimming)
process trimmomatic {
    container 'biocontainers/trimmomatic:v0.39'
    input:
    path reads

    output:
    path 'zygmukin_trimmed_reads'

    script:
    """
    trimmomatic PE -threads 4 ${reads[0]} ${reads[1]} \
        zygmukin_trimmed_reads_1P.fastq.gz zygmukin_trimmed_reads_1U.fastq.gz \
        zygmukin_trimmed_reads_2P.fastq.gz zygmukin_trimmed_reads_2U.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Process 3: Trinity Assembly (De novo assembly)
process trinity_assembly {
    container 'trinityrnaseq/trinityrnaseq:2.13.2'
    input:
    path reads

    output:
    path 'zygmukin_trinity_output'

    script:
    """
    Trinity --seqType fq --max_memory 50G --left ${reads[0]} --right ${reads[1]} --CPU 4 --output zygmukin_trinity_output
    """
}

// Process 4: RSEM (RNA-Seq by Expectation-Maximization)
process rsem_analysis {
    container 'biocontainers/rsem:v1.3.1'
    input:
    path 'zygmukin_trinity_output/Trinity.fasta'   // Transcriptome for RSEM reference
    path 'zygmukin_trimmed_reads_1P.fastq.gz'       // Trimmed paired-end reads
    path 'zygmukin_trimmed_reads_2P.fastq.gz'       // Trimmed paired-end reads

    output:
    path 'zygmukin_rsem_results'

    script:
    """
    rsem-calculate-expression --paired-end \
        --bowtie2 --alignments \
        --output-genome-bam --estimate-rspd \
        --strand-specific rsem_output \
        ${reads[0]} ${reads[1]} zygmukin_trinity_output/Trinity.fasta
    """
}

// Process 5: TransDecoder (ORF prediction)
process transdecoder_prediction {
    container 'biocontainers/transdecoder:v5.5.0'
    input:
    path 'zygmukin_trinity_output/Trinity.fasta'

    output:
    path 'zygmukin_transdecoder_output'

    script:
    """
    TransDecoder.LongOrfs -t Trinity.fasta
    TransDecoder.Predict -t Trinity.fasta
    mv Trinity.fasta.transdecoder.* zygmukin_transdecoder_output
    """
}

// Process 6: Trinotate Annotation
process trinnotate_annotation {
    container 'trinotate/trinotate:latest'
    input:
    path 'zygmukin_trinity_output/Trinity.fasta'
    path 'counts.txt'

    output:
    path 'zygmukin_trinotate_results'

    script:
    """
    ./trinotate_annotation.sh zygmukin_trinity_output/Trinity.fasta counts.txt zygmukin_trinotate_results
    """
}

// Process 7: KOBAS Annotation (Functional Annotation)
process kobas_annotation {
    container 'agbase/kobas:3.0.3_3'
    input:
    path 'zygmukin_transdecoder_output/Trinity.fasta.transdecoder.pep'

    output:
    path 'zygmukin_kobas_results'

    script:
    """
    sudo docker run --rm -v $(pwd):/work-dir agbase/kobas:3.0.3_3 \
        -a -i /work-dir/zygmukin_transdecoder_output/Trinity.fasta.transdecoder.pep \
        -s ath -t fasta:nuc -o /work-dir/zygmukin_kobas_results
    """
}

// Process 8: DESeq2 Analysis (Differential Expression)
process deseq2_analysis {
    container 'biocontainers/r-ver:v3.6.1'
    input:
    path 'counts.txt'

    output:
    path 'zygmukin_deseq2_results.csv'

    script:
    """
    Rscript deseq2_analysis.R counts.txt zygmukin_deseq2_results.csv
    """
}

// Process 9: Limma-voom Analysis (Differential Expression)
process limma_voom_analysis {
    container 'biocontainers/r-ver:v3.6.1'
    input:
    path 'counts.txt'

    output:
    path 'zygmukin_limma_voom_results.csv'

    script:
    """
    Rscript limma_voom_analysis.R counts.txt zygmukin_limma_voom_results.csv
    """
}

// Process 10: edgeR Analysis (Differential Expression)
process edger_analysis {
    container 'biocontainers/r-ver:v3.6.1'
    input:
    path 'counts.txt'

    output:
    path 'zygmukin_edger_results.csv'

    script:
    """
    Rscript edger_analysis.R counts.txt zygmukin_edger_results.csv
    """
}

workflow {
    fastqc(params.reads)  // Quality control check
    trimmomatic(params.reads) | trinity_assembly()  // Read trimming → Transcriptome assembly
    trinity_assembly() | rsem_analysis()  // RNA-Seq quantification with RSEM
    trinity_assembly() | transdecoder_prediction()  // Find coding sequences
    trinity_assembly() | trinnotate_annotation()  // Transcript annotation
    transdecoder_prediction() | kobas_annotation()  // Functional annotation of ORFs
    deseq2_analysis() | limma_voom_analysis()  // Differential expression analysis with limma-voom
    deseq2_analysis() | edger_analysis()  // Differential expression analysis with edgeR
}

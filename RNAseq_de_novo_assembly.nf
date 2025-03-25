// Define input parameters and set the appropriate paths for your files
params.reads = 'C:/tools_phd/zygmukin_reads/*.fastq.gz'
params.outdir = 'C:/tools_phd/zygmukin_results'
params.transcriptome = 'C:/tools_phd/zygmukin_trinity_output/Trinity.fasta'  // Path to Trinity assembly

// Process 1: FastQC (Quality Control)
process fastqc {
    container 'biocontainers/fastqc:v0.11.9'
    input:
    path reads

    output:
    path 'C:/tools_phd/zygmukin_fastqc_reports'

    script:
    """
    fastqc -o C:/tools_phd/zygmukin_fastqc_reports $reads
    """
}

// Process 2: Trimmomatic (Read Trimming)
process trimmomatic {
    container 'biocontainers/trimmomatic:v0.39'
    input:
    path reads

    output:
    path 'C:/tools_phd/zygmukin_trimmed_reads'

    script:
    """
    trimmomatic PE -threads 4 ${reads[0]} ${reads[1]} \
        C:/tools_phd/zygmukin_trimmed_reads_1P.fastq.gz C:/tools_phd/zygmukin_trimmed_reads_1U.fastq.gz \
        C:/tools_phd/zygmukin_trimmed_reads_2P.fastq.gz C:/tools_phd/zygmukin_trimmed_reads_2U.fastq.gz \
        ILLUMINACLIP:C:/tools_phd/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Process 3: Trinity Assembly (De novo assembly)
process trinity_assembly {
    container 'trinityrnaseq/trinityrnaseq:2.13.2'
    input:
    path reads

    output:
    path 'C:/tools_phd/zygmukin_trinity_output'

    script:
    """
    Trinity --seqType fq --max_memory 50G --left ${reads[0]} --right ${reads[1]} --CPU 4 --output C:/tools_phd/zygmukin_trinity_output
    """
}

// Process 4: RSEM (RNA-Seq by Expectation-Maximization)
process rsem_analysis {
    container 'biocontainers/rsem:v1.3.1'
    input:
    path 'C:/tools_phd/zygmukin_trinity_output/Trinity.fasta'   // Transcriptome for RSEM reference
    path 'C:/tools_phd/zygmukin_trimmed_reads_1P.fastq.gz'       // Trimmed paired-end reads
    path 'C:/tools_phd/zygmukin_trimmed_reads_2P.fastq.gz'       // Trimmed paired-end reads

    output:
    path 'C:/tools_phd/zygmukin_rsem_results'

    script:
    """
    rsem-calculate-expression --paired-end \
        --bowtie2 --alignments \
        --output-genome-bam --estimate-rspd \
        --strand-specific C:/tools_phd/rsem_output \
        ${reads[0]} ${reads[1]} C:/tools_phd/zygmukin_trinity_output/Trinity.fasta
    """
}

// Process 5: TransDecoder (ORF prediction)
process transdecoder_prediction {
    container 'biocontainers/transdecoder:v5.5.0'
    input:
    path 'C:/tools_phd/zygmukin_trinity_output/Trinity.fasta'

    output:
    path 'C:/tools_phd/zygmukin_transdecoder_output'

    script:
    """
    TransDecoder.LongOrfs -t C:/tools_phd/Trinity.fasta
    TransDecoder.Predict -t C:/tools_phd/Trinity.fasta
    mv C:/tools_phd/Trinity.fasta.transdecoder.* C:/tools_phd/zygmukin_transdecoder_output
    """
}

// Process 6: Trinotate Annotation
process trinnotate_annotation {
    container 'trinotate/trinotate:latest'
    input:
    path 'C:/tools_phd/zygmukin_trinity_output/Trinity.fasta'
    path 'C:/tools_phd/counts.txt'

    output:
    path 'C:/tools_phd/zygmukin_trinotate_results'

    script:
    """
    ./trinotate_annotation.sh C:/tools_phd/zygmukin_trinity_output/Trinity.fasta C:/tools_phd/counts.txt C:/tools_phd/zygmukin_trinotate_results
    """
}

// Process 7: KOBAS Annotation (Functional Annotation)
process kobas_annotation {
    container 'agbase/kobas:3.0.3_3'
    input:
    path 'C:/tools_phd/zygmukin_transdecoder_output/Trinity.fasta.transdecoder.pep'

    output:
    path 'C:/tools_phd/zygmukin_kobas_results'

    script:
    """
    sudo docker run --rm -v C:/tools_phd:/work-dir agbase/kobas:3.0.3_3 \
        -a -i /work-dir/zygmukin_transdecoder_output/Trinity.fasta.transdecoder.pep \
        -s ath -t fasta:nuc -o /work-dir/zygmukin_kobas_results
    """
}

// Process 8: DESeq2 Analysis (Differential Expression)
process deseq2_analysis {
    container 'biocontainers/r-ver:v3.6.1'
    input:
    path 'C:/tools_phd/counts.txt'

    output:
    path 'C:/tools_phd/zygmukin_deseq2_results.csv'

    script:
    """
    Rscript C:/tools_phd/deseq2_analysis.R C:/tools_phd/counts.txt C:/tools_phd/zygmukin_deseq2_results.csv
    """
}

// Process 9: Limma-voom Analysis (Differential Expression)
process limma_voom_analysis {
    container 'biocontainers/r-ver:v3.6.1'
    input:
    path 'C:/tools_phd/counts.txt'

    output:
    path 'C:/tools_phd/zygmukin_limma_voom_results.csv'

    script:
    """
    Rscript C:/tools_phd/limma_voom_analysis.R C:/tools_phd/counts.txt C:/tools_phd/zygmukin_limma_voom_results.csv
    """
}

// Process 10: edgeR Analysis (Differential Expression)
process edger_analysis {
    container 'biocontainers/r-ver:v3.6.1'
    input:
    path 'C:/tools_phd/counts.txt'

    output:
    path 'C:/tools_phd/zygmukin_edger_results.csv'

    script:
    """
    Rscript C:/tools_phd/edger_analysis.R C:/tools_phd/counts.txt C:/tools_phd/zygmukin_edger_results.csv
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

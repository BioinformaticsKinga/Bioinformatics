RNA-Seq Analysis Pipeline using Nextflow


Overview

This repository contains a Nextflow pipeline for RNA-Seq analysis, including quality control, read trimming, transcriptome assembly, expression quantification, and differential expression analysis.

Pipeline Steps

The pipeline consists of the following processes:

FastQC - Quality control of raw reads

Trimmomatic - Read trimming to remove adapters and low-quality sequences

Trinity - De novo transcriptome assembly

RSEM - Expression quantification

TransDecoder - Open Reading Frame (ORF) prediction

Trinotate - Transcriptome annotation

KOBAS - Functional annotation of predicted proteins

DESeq2, Limma-voom, edgeR - Differential expression analysis

Transcriptome Assembly

The de novo transcriptome assembly in this pipeline is performed using Trinity, a tool specifically designed for reconstructing transcriptomes from RNA-Seq data in the absence of a reference genome. The assembly process involves three main stages:

Inchworm: Assembles contigs based on k-mer frequency.

Chrysalis: Clusters contigs into potential transcript families and constructs de Bruijn graphs.

Butterfly: Resolves alternative splicing and transcript variations to generate final transcript sequences.

This step is critical for downstream analyses, including differential expression and functional annotation.

For implementation details, refer to the Nextflow script: RNAseq_de_novo_assembly.nf.

Installation

Prerequisites

Nextflow installed

Docker or Singularity installed (for containerized execution)

Reference files (e.g., Trinity.fasta, adapter sequences for trimming)

Clone the Repository

git clone https://github.com/your-username/your-repo.git
cd your-repo

Running the Pipeline

With Default Parameters

nextflow run main.nf -profile docker

With Custom Parameters

Modify nextflow.config or provide parameters via command line:

nextflow run main.nf --reads '/path/to/reads/*.fastq.gz' --outdir '/path/to/output' --transcriptome '/path/to/Trinity.fasta'

Configuration

The pipeline can be configured via nextflow.config. Ensure that paths to necessary reference files and tools are correctly set.

Output Files

The pipeline generates several output directories:

zygmukin_fastqc_reports/ - FastQC reports

zygmukin_trimmed_reads/ - Trimmed reads

zygmukin_trinity_output/ - Trinity transcriptome assembly

zygmukin_rsem_results/ - Expression quantification results

zygmukin_transdecoder_output/ - ORF predictions

zygmukin_trinotate_results/ - Functional annotation

zygmukin_kobas_results/ - Pathway annotations

zygmukin_deseq2_results.csv - DESeq2 differential expression results

zygmukin_limma_voom_results.csv - Limma-voom results

zygmukin_edger_results.csv - edgeR results

Troubleshooting

Ensure all required reference files exist.

Check Docker/Singularity setup if using containers.

Run with -resume if you want to continue from a failed step.

License

This project is licensed under the MIT License.

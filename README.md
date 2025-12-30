# RNA-seq Analysis Pipeline: DUS1L Knockout Study

End-to-end RNA-seq differential expression pipeline using Nextflow, Docker, and AWS.

## Project Overview

This pipeline analyzes RNA-seq data comparing **wild-type (WT) vs DUS1L knockout (KO)** cells to identify differentially expressed genes.

**Dataset**: [GSE312810](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312810)
- 6 samples (3 WT, 3 KO)
- Paired-end Illumina sequencing
- Human genome (GRCh38)

## Pipeline Steps
```
FASTQ Files
    ↓
FastQC (Quality Control)
    ↓
Fastp (Adapter Trimming & Filtering)
    ↓
Salmon (Pseudo-alignment & Quantification)
    ↓
DESeq2 (Differential Expression Analysis)
    ↓
Visualizations (PCA, Volcano, Heatmaps)
```

## Key Results

- **22,213 genes** analyzed after filtering
- **5,155 significantly DE genes** (FDR < 0.05)
  - 2,533 upregulated in KO
  - 2,622 downregulated in KO

## Technologies Used

- **Workflow Manager**: Nextflow (DSL2)
- **Containerization**: Docker
- **Cloud Platform**: AWS EC2 (r6i.2xlarge)
- **Tools**: FastQC, Fastp, Salmon, DESeq2
- **Languages**: Groovy (Nextflow), R (DESeq2)

## Repository Structure
```
rna-seq-project/
├── workflow/
│   ├── main.nf              # Nextflow pipeline
│   └── config/
│       └── samplesheet.csv  # Sample metadata
├── containers/
│   └── Dockerfile           # Docker image with all tools
├── scripts/
│   └── deseq2_analysis.R    # Differential expression analysis
├── config/
├── results/
│   ├── fastqc/              # QC reports
│   ├── fastp/               # Trimming reports
│   ├── salmon/              # Quantification results
│   └── deseq2/              # DE analysis & plots
└── README.md
```

## Quick Start

### Prerequisites
- Nextflow ≥ 25.10
- Docker
- AWS account (optional, for cloud execution)

### Local Execution
```bash
# Clone repository
git clone https://github.com/ajay-shenoy4/RNA-Seq-Pipeline.git
cd RNA-Seq-Pipeline

# Build Docker image
docker build -t rna-seq-pipeline:v1 containers/

# Download reference genome (example)
cd data/reference
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz
gunzip gencode.v46.transcripts.fa.gz

# Run pipeline
nextflow run workflow/main.nf -c nextflow.config
```

### AWS EC2 Execution
```bash
# Launch EC2 instance
aws ec2 run-instances \
  --image-id ami-0c02fb55c47d7f853 \
  --instance-type r6i.2xlarge \
  --key-name your-key \
  --security-group-ids sg-xxx

# SSH and run pipeline
ssh -i "your-key.pem" ubuntu@YOUR-IP
# ... (setup and run as above)
```

**Cost**: ~$0.50/hour (r6i.2xlarge)
**Runtime**: ~45 minutes for full pipeline

## Key Visualizations

### PCA Plot
Shows clear separation between WT and KO samples along PC1.

### Volcano Plot
Identifies 5,155 DE genes with strong effect sizes (log2FC > 1).

### MA Plot
Distribution of differential expression across expression levels.

### Heatmap
Top 50 most significantly DE genes cluster by condition.

## Output Files

### Quantification
- `results/salmon/*/quant.sf` - Gene-level abundance estimates

### Differential Expression
- `results/deseq2/salmon_all_genes.csv` - Full DE results
- `results/deseq2/salmon_significant_genes.csv` - FDR < 0.05

### Visualizations
- `results/deseq2/salmon_pca_plot.pdf`
- `results/deseq2/salmon_volcano_plot.pdf`
- `results/deseq2/salmon_ma_plot.pdf`
- `results/deseq2/salmon_top50_genes_heatmap.pdf`

## Configuration

Edit `nextflow.config` to customize:
- CPU/memory allocation
- Enable STAR alignment (`run_star = true`)
- Change reference genome paths

## Citation

Dataset: GSE312810 - DUS1L knockout RNA-seq study

## License

MIT License

## Author

Ajay Shenoy
- GitHub: [@ajay-shenoy4](https://github.com/ajay-shenoy4)
- Project: [RNA-Seq-Pipeline](https://github.com/ajay-shenoy4/RNA-Seq-Pipeline)

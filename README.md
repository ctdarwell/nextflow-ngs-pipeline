## Overview

This pipeline performs:
- Quality control (FastQC)
- Read alignment (BWA)
- BAM conversion (samtools)
- Alignment statistics (flagstat)
- Aggregated reporting (MultiQC)

NB as is, for toy dataset (config file allows scalable resources: `mem`, `ncpus` etc)

## Expected output

After running, you should see:

- `results/fastqc/` → per-sample QC reports  
- `results/bam/` → aligned BAM files  
- `results/flagstat/` → mapping statistics  
- `results/multiqc/multiqc_report.html` → aggregated report  

# Usage
```bash
nextflow run main.nf
```

# Inputs
- FASTQ files: data/Hd4*_{1,2}.fastq.gz
- Reference genome: data/*.fa

# Outputs
- results/fastqc/
- results/bam/
- results/flagstat/
- results/multiqc/

# Requirements
- Nextflow
- Docker

## Initial setup: install dependencies

```bash
# Install Java, Nextflow, Docker
sudo apt update
sudo apt install openjdk-17-jdk -y

java -version

curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/
export PATH=$HOME/bin:$PATH

sudo apt install docker.io -y
sudo usermod -aG docker $USER

newgrp docker
docker run hello-world

sudo apt install bwa

```

## STEP 2: ADD DATA

```bash
# 1. Create a directory for the test data
mkdir -p data && cd data

# 2. Download the Human Reference Genome (GRCh38)
# We use the 'primary assembly' which is the standard for mapping
echo "Downloading Human Reference Genome (GRCh38)..."
wget -c ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# 3. Download 4 Human Samples (from 1000 Genomes Project via ENA)
# These are paired-end (.1 and .2) files.
# Using small accessions for demonstration; these are real human data.
samples=("ERR008901" "ERR008902" "ERR008903" "ERR008904")

echo "Downloading 4 human FASTQ samples..."
for id in "${samples[@]}"; do
    echo "Getting sample $id..."
    # Format the ENA FTP path based on accession subfolders
    # Example: /vol1/fastq/ERR008/ERR008901/ERR008901_1.fastq.gz
    prefix=${id:0:6}
    wget -c "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$prefix/$id/${id}_1.fastq.gz"
    wget -c "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/$prefix/$id/${id}_2.fastq.gz"
done

echo "Download complete."

# 4. Reduce to toy data size (100000 reads)
for id in $(ls ERR*gz); do
    echo "$id"
	zcat "$id" | head -n 400000 | gzip > "Hd400k_$id"

done


# 5. unzip reference genome
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# 6. create index files for reference
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa

```


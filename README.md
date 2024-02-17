# Variant calling pipeline for Pacbio Hifi reads

## Prerequisites
Ensure that the following tools are installed:

- bcftools
- Docker
- pbmm2
- pbsv
- samtools
- seqkit
- whatshap

> [!IMPORTANT]
> Please note that Docker is required to run certain tools like deepvariant and ensembl-vep which are executed within Docker containers.

## Running the Pipeline

The pipeline requires the following inputs:
- `reads_fastq_gz`
- `file_label`-
- `genome_ref`
- `vep_cache`

> [!WARNING]
> If you intend to run the entire pipeline, including the annotation step, please specify the VEP cache directory. If the annotation step is not required, you can run the appropriate target as shown in the help message.

To run the pipeline, use the following command:

```bash
make reads_fastq_gz=./test_data/HG002-downsample-0.1.fastq.gz file_label=HG002 genome_ref=./test_data/Homo_sapiens.GRCh38.dna.chromosome.22.fa vep_cache=./vep_cache threads=4
```
## Complete usage

Type the following command to see the full help message.

```{bash}
make
```

Output:
```{bash}
Usage: make [TARGET] reads_fastq_gz=<path> file_label=<label> genome_ref=<path> vep_cache=<path> threads=<integer>

Targets: (default: all)
  alignReads callSmallVars callStructuralVars phaseVars annotateVars

Options:
  reads_fastq_gz=<path>    Path to the reads_fastq_gz file (required)
  file_label=<label>       Label for the file (required)
  genome_ref=<path>        Path to the genome_ref file (required)
  vep_cache=<path>         Path to the vep_cache file (required)
  output=<path>            Path to the output directory (optional). Default: /home/anand/Documents/aspire-files/data-documents/make-pipeline
  threads=<integer>        Number of threads (optional). Default: 12
  -n                       dry-run [simulation of the build process without actually executing the commands] (optional)

Example:

  >To run all the steps:
  make reads_fastq_gz=./test_data/HG002-downsample-0.1.fastq.gz file_label=HG002 genome_ref=./test_data/Homo_sapiens.GRCh38.dna.chromosome.22.fa vep_cache=./vep_cache threads=4

  >Dry run:
  make reads_fastq_gz=./test_data/HG002-downsample-0.1.fastq.gz file_label=HG002 genome_ref=./test_data/Homo_sapiens.GRCh38.dna.chromosome.22.fa vep_cache=./vep_cache threads=4 -n

  >To run a single target:
  make alignReads reads_fastq_gz=./test_data/HG002-downsample-0.1.fastq.gz file_label=HG002 genome_ref=./test_data/Homo_sapiens.GRCh38.dna.chromosome.22.fa vep_cache=./vep_cache threads=4
```

## Example

### Downloading the test data

Before running the Makefile, you need to download the necessary data. This can be done using the `download.sh` script provided in the `test_data` directory. This script will download a fastq file and a reference genome if they do not already exist in the current directory. It will also downsample the fastq file.

To run the script, navigate to the directory containing the script and use the following command:

```{bash}
bash download.sh
```

After running the script, you should see the following files in your current directory:

- `HG002.fastq.gz`: The original fastq file.
- `Homo_sapiens.GRCh38.dna.primary_assembly.fa`: The reference genome.
- `HG002-downsample-0.1.fastq.gz`: The downsampled fastq file.

### Execution

You can then use the files downloaded in the previous step as inputs to the `Makefile`. For example:

```{bash}
cd make-pipeline

make reads_fastq_gz=test_data/HG002-downsample-0.1.fastq.gz genome_ref=test_data/Homo_sapiens.GRCh38.dna.primary_assembly.fa file_label=demo
```

This command will run the Makefile with the downsampled fastq file and the reference genome, and it will label the output files with `results-demo`.

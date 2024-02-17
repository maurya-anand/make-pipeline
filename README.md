# make-pipeline

## Prerequisites
Ensure that the following tools are installed:

- bcftools
- Docker
- pbmm2
- pbsv
- samtools
- seqkit
- whatshap

> [!TIP]
> Please note that Docker is required to run certain tools like deepvariant and ensembl-vep which are executed within Docker containers.


## Example usage

### Downloading the Data

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

### Complete usage

Type the following command to see the full help message.

```{bash}
make
```

Output:
```{bash}
Usage: make [TARGET] reads_fastq_gz=<path> file_label=<label> [other options]

Targets: (default: all)
  alignReads callSmallVars callStructuralVars phaseVars annotateVars

Options:
  reads_fastq_gz=<path>    Path to the reads_fastq_gz file (required)
  file_label=<label>       Label for the file (required)
  genome_ref=<path>        Path to the genome_ref file (optional)
  vep_cache=<path>         Path to the vep_cache file (optional)
  -n                       dry-run [simulation of the build process without actually executing the commands] (optional)

Example:
  make reads_fastq_gz=g2.fastq.gz file_label=test -n
  make reads_fastq_gz=g2.fastq.gz file_label=test -n
  make alignReads reads_fastq_gz=g2.fastq.gz file_label=test -n
```
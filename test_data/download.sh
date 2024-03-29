cd "$(dirname "$0")"
SCRIPT_DIR="$(pwd)"

if [ ! -f "HG002.fastq.gz" ]; then
    curl https://downloads.pacbcloud.com/public/dataset/HiFiTE_SqIIe/Oct_2022/TwistAllianceLongReadPGx/fastq/HG002.fastq.gz --output HG002.fastq.gz
fi
if [ ! -f "Homo_sapiens.GRCh38.dna.chromosome.22.fa" ]; then
    curl https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz --output Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz; \
    gunzip Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
fi
if [ ! -f "HG002-downsample-0.1.fastq.gz" ]; then
    # downsampling the downloaded fastq
    seqtk sample HG002.fastq.gz 0.1 > HG002-downsample-0.1.fastq
    gzip HG002-downsample-0.1.fastq
fi
printf "\nTest data downloaded in: $SCRIPT_DIR\n";
printf "\nIf you intend to run the entire pipeline, including the annotation step, please specify the VEP cache directory. If the annotation step is not required, you can run the appropriate target as shown in the help message.\n";
printf "\nTry:\nmake reads_fastq_gz=$SCRIPT_DIR/test_data/HG002-downsample-0.1.fastq.gz file_label=HG002 genome_ref=$SCRIPT_DIR/test_data/Homo_sapiens.GRCh38.dna.chromosome.22.fa vep_cache=<path_to_vep_cache> threads=4\n";
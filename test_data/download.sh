cd "$(dirname "$0")"
SCRIPT_DIR="$(pwd)"

if [ ! -f "HG002.fastq.gz" ]; then
    curl https://downloads.pacbcloud.com/public/dataset/HiFiTE_SqIIe/Oct_2022/TwistAllianceLongReadPGx/fastq/HG002.fastq.gz --output HG002.fastq.gz
fi
if [ ! -f "Homo_sapiens.GRCh38.dna.chromosome.22.fa" ]; then
    curl https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz --output Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz; \
    gunzip Homo_sapiens.GRCh38.dna.chromosome.22.fa.gz
fi

# downsampling the downloaded fastq
seqtk sample HG002.fastq.gz 0.1 > HG002-downsample-0.1.fastq
gzip HG002-downsample-0.1.fastq

echo "Test data downloaded in: $SCRIPT_DIR"

echo "Try: make reads_fastq_gz=test_data/HG002-downsample-0.1.fastq.gz genome_ref=test_data/Homo_sapiens.GRCh38.dna.chromosome.22.fa file_label=demo"
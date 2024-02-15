cd "$(dirname "$0")"
SCRIPT_DIR="$(pwd)"

if [ ! -f "HG002.fastq.gz" ]; then
    curl https://downloads.pacbcloud.com/public/dataset/HiFiTE_SqIIe/Oct_2022/TwistAllianceLongReadPGx/fastq/HG002.fastq.gz --output HG002.fastq.gz
fi
if [ ! -f "Homo_sapiens.GRCh38.dna.primary_assembly.fa" ]; then
    curl https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz --output Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz; \
    gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
fi

echo "Test data downloaded in: $SCRIPT_DIR"
# Mandatory inputs
reads_fastq_gz ?=
file_label ?=
genome_ref ?=
vep_cache ?=

# Optional inputs
output ?= $(PWD)/results-${file_label}
threads ?= 4

ifeq ($(reads_fastq_gz),)
.DEFAULT_GOAL := help
endif

ifeq ($(file_label),)
.DEFAULT_GOAL := help
endif

ifeq ($(genome_ref),)
.DEFAULT_GOAL := help
endif

ifneq ($(reads_fastq_gz),)
ifneq ($(file_label),)
ifneq ($(genome_ref),)
.DEFAULT_GOAL := all
endif
endif
endif

STEP1_DIR := $(output)/01_alignReads
STEP2_DIR := $(output)/02_callSmallVars
STEP3_DIR := $(output)/03_callStructuralVars
STEP4_DIR := $(output)/04_phaseVars
STEP5_DIR := $(output)/05_annotateVars
TEMP_DIR := $(output)/temp

all: alignReads callSmallVars callStructuralVars phaseVars annotateVars
	@printf "\nPipeline completed.\n\nResults dir: ${output}\n";

help:
	@echo "Usage: make [TARGET] reads_fastq_gz=<path> file_label=<label> genome_ref=<path> vep_cache=<path> threads=<integer>"
	@echo ""
	@echo "Targets: (default: all)"
	@echo "  alignReads callSmallVars callStructuralVars phaseVars annotateVars"
	@echo ""
	@echo "Options:"
	@echo "  reads_fastq_gz=<path>    Path to the reads_fastq_gz file (required)"
	@echo "  file_label=<label>       Label for the file (required)"
	@echo "  genome_ref=<path>        Path to the genome_ref file (required)"
	@echo "  vep_cache=<path>         Path to the vep_cache file (required)"
	@echo "  output=<path>            Path to the output directory (optional). Default: $(PWD)"
	@echo "  threads=<integer>        Number of threads (optional). Default: 12"
	@echo "  -n                       dry-run [simulation of the build process without actually executing the commands] (optional)"
	@echo ""
	@echo "Example:"
	@echo ""
	@echo "  >To run all the steps:"
	@echo "  make reads_fastq_gz=./test_data/HG002-downsample-0.1.fastq.gz file_label=HG002 genome_ref=./test_data/Homo_sapiens.GRCh38.dna.chromosome.22.fa vep_cache=./vep_cache threads=4"
	@echo ""
	@echo "  >Dry run:"
	@echo "  make reads_fastq_gz=./test_data/HG002-downsample-0.1.fastq.gz file_label=HG002 genome_ref=./test_data/Homo_sapiens.GRCh38.dna.chromosome.22.fa vep_cache=./vep_cache threads=4 -n"
	@echo ""
	@echo "  >To run a single target:"
	@echo "  make alignReads reads_fastq_gz=./test_data/HG002-downsample-0.1.fastq.gz file_label=HG002 genome_ref=./test_data/Homo_sapiens.GRCh38.dna.chromosome.22.fa vep_cache=./vep_cache threads=4"

alignReads:
	@if [ -z "$(reads_fastq_gz)" ]; then \
		echo "reads_fastq_gz is not set. Run 'make help' for usage instructions."; \
		exit 1; \
	fi
	@if [ -z "$(file_label)" ]; then \
		echo "file_label is not set. Run 'make help' for usage instructions."; \
		exit 1; \
	fi
	@if [ -z "$(genome_ref)" ]; then \
		echo "genome_ref is not set. Run 'make help' for usage instructions."; \
		exit 1; \
	fi
	@if [ ! -f ${STEP1_DIR}/${file_label}_aligned.bam ]; then \
		mkdir -p $(STEP1_DIR) $(TEMP_DIR); \
		chmod -R a+w $(STEP1_DIR); \
		cp ${genome_ref} $(TEMP_DIR)/genome_reference.fasta; \
		samtools faidx $(TEMP_DIR)/genome_reference.fasta -o $(TEMP_DIR)/genome_reference.fasta.fai; \
		gunzip -c ${reads_fastq_gz} > $(TEMP_DIR)/${file_label}.hifi_reads.fastq; \
		pbmm2 align \
			--preset hifi \
			--sort \
			$(TEMP_DIR)/genome_reference.fasta \
			$(TEMP_DIR)/${file_label}.hifi_reads.fastq \
			${STEP1_DIR}/${file_label}_aligned.bam \
			--rg "@RG\tID:${file_label}\tSM:${file_label}" \
			--log-level DEBUG \
			--log-file ${STEP1_DIR}/${file_label}_alignment.log; \
		seqkit stats -a -T $(TEMP_DIR)/${file_label}.hifi_reads.fastq > ${STEP1_DIR}/${file_label}_fastq_stats.tab; \
		samtools depth ${STEP1_DIR}/${file_label}_aligned.bam -o ${STEP1_DIR}/${file_label}_depth.txt; \
    fi
	@printf "\nFinished aligning reads to reference.\n";

callSmallVars: alignReads
	@if [ ! -f ${STEP2_DIR}/${file_label}_ALL_variants.vcf.gz ]; then \
		mkdir -p $(STEP2_DIR); \
		chmod -R a+w $(STEP2_DIR); \
		if [ $$(samtools view ${STEP1_DIR}/${file_label}_aligned.bam | wc -l) -eq 0 ]; then \
			touch ${STEP2_DIR}/${file_label}_ALL_variants.vcf.gz ${STEP2_DIR}/${file_label}_all_variants.visual_report.html ${STEP2_DIR}/${file_label}_PASS_variants.vcf.gz ${STEP2_DIR}/${file_label}_PASS_NORM_variants.vcf.gz; \
		else \
			docker run \
				-v $(STEP1_DIR):/input \
				-v $(TEMP_DIR):/temp \
				-v $(STEP2_DIR):/output \
				google/deepvariant:1.5.0 \
				/opt/deepvariant/bin/run_deepvariant \
				--model_type PACBIO \
				--num_shards ${threads} \
				--ref /temp/genome_reference.fasta \
				--reads /input/${file_label}_aligned.bam \
				--output_vcf /output/${file_label}_ALL_variants.vcf.gz >> ${STEP2_DIR}/${file_label}_variant_calling.log 2>&1; \
			bcftools view -f PASS ${STEP2_DIR}/${file_label}_ALL_variants.vcf.gz -Oz -o ${STEP2_DIR}/${file_label}_PASS_variants.vcf.gz >> ${STEP2_DIR}/${file_label}_variant_calling.log 2>&1; \
			bcftools norm ${STEP2_DIR}/${file_label}_PASS_variants.vcf.gz -f ${TEMP_DIR}/genome_reference.fasta -m -any -Oz -o ${STEP2_DIR}/${file_label}_PASS_NORM_variants.vcf.gz >> ${STEP2_DIR}/${file_label}_variant_calling.log 2>&1; \
			tabix -f -p vcf ${STEP2_DIR}/${file_label}_PASS_NORM_variants.vcf.gz; \
		fi; \
	fi
	@printf "\nFinished calling small variants.\n";

callStructuralVars: alignReads
	@if [ ! -f ${STEP3_DIR}/${file_label}_PASS_NORM_structural_variants.vcf ]; then \
		mkdir -p $(STEP3_DIR); \
		chmod -R a+w $(STEP3_DIR); \
		pbsv discover \
			${STEP1_DIR}/${file_label}_aligned.bam \
			${STEP3_DIR}/${file_label}_aligned.svsig.gz >> ${STEP3_DIR}/${file_label}_structural_variant_calling.log 2>&1; \
		pbsv call \
			${TEMP_DIR}/genome_reference.fasta \
			${STEP3_DIR}/${file_label}_aligned.svsig.gz \
			${STEP3_DIR}/${file_label}_ALL_structural_variants.vcf >> ${STEP3_DIR}/${file_label}_structural_variant_calling.log 2>&1; \
		bcftools view -f PASS ${STEP3_DIR}/${file_label}_ALL_structural_variants.vcf -Ov -o ${STEP3_DIR}/${file_label}_PASS_structural_variants.vcf >> ${STEP3_DIR}/${file_label}_structural_variant_calling.log 2>&1; \
		bcftools norm ${STEP3_DIR}/${file_label}_PASS_structural_variants.vcf -f ${TEMP_DIR}/genome_reference.fasta -m -any -Ov -o ${STEP3_DIR}/${file_label}_PASS_NORM_structural_variants.vcf >> ${STEP3_DIR}/${file_label}_structural_variant_calling.log 2>&1; \
	fi
	@printf "\nFinished calling structural variants.\n";

phaseVars: callSmallVars
	@if [ ! -f ${STEP4_DIR}/${file_label}_PASS_NORM_PHASED_variants.vcf ]; then \
		mkdir -p $(STEP4_DIR); \
		chmod -R a+w $(STEP4_DIR); \
		if [ $$(zcat ${STEP2_DIR}/${file_label}_PASS_NORM_variants.vcf.gz | grep -v "#" | wc -l) -eq 0 ]; then \
			touch ${STEP4_DIR}/${file_label}_PASS_NORM_PHASED_variants.vcf ${STEP4_DIR}/${file_label}_variant_phasing.log; \
		else \
			whatshap phase \
				${STEP2_DIR}/${file_label}_PASS_NORM_variants.vcf.gz \
				${STEP1_DIR}/${file_label}_aligned.bam \
				--reference=${TEMP_DIR}/genome_reference.fasta \
				--indels \
				--output ${STEP4_DIR}/${file_label}_PASS_NORM_PHASED_variants.vcf >> ${STEP4_DIR}/${file_label}_variant_phasing.log 2>&1; \
			whatshap stats \
				${STEP4_DIR}/${file_label}_PASS_NORM_PHASED_variants.vcf \
				--tsv=${STEP4_DIR}/${file_label}_variant_phasing.log >> ${STEP4_DIR}/${file_label}_variant_phasing.log 2>&1; \
		fi; \
	fi
	@printf "\nFinished phasing small variants.\n";

annotateVars: phaseVars
	@if [ -z "$(vep_cache)" ]; then \
		echo "vep_cache is not set. Run 'make help' for usage instructions."; \
		exit 1; \
	fi
	@if [ ! -f ${STEP5_DIR}/${file_label}_PASS_NORM_PHASED_ANNOTATED_variants.vcf ]; then \
		mkdir -p $(STEP5_DIR); \
		chmod -R a+w $(STEP5_DIR); \
		if [ $$(grep -v "#" ${STEP4_DIR}/${file_label}_PASS_NORM_PHASED_variants.vcf | wc -l) -eq 0 ]; then \
			touch ${STEP5_DIR}/${file_label}_PASS_NORM_PHASED_ANNOTATED_variants.vcf ${STEP5_DIR}/${file_label}_vep_stats.log; \
		else \
			docker run \
				-v ${STEP4_DIR}:/input \
				-v $(TEMP_DIR):/temp \
				-v ${vep_cache}:/vep_cache \
				-v $(STEP5_DIR):/output \
				ensemblorg/ensembl-vep:release_110.1 \
				perl /opt/vep/src/ensembl-vep/vep --force_overwrite \
				--input_file /input/${file_label}_PASS_NORM_PHASED_variants.vcf \
				--vcf \
				--output_file /output/${file_label}_PASS_NORM_PHASED_ANNOTATED_variants.vcf \
				--stats_file /output/${file_label}_vep_stats.log \
				--stats_text \
				--cache \
				--dir_cache /vep_cache \
				--fork ${threads} \
				--fasta /temp/genome_reference.fasta \
				--numbers --offline --hgvs --shift_hgvs 0 --terms SO --symbol \
				--sift b --polyphen b --total_length --ccds --canonical --biotype \
				--protein --xref_refseq --mane --pubmed --af --max_af --af_1kg --af_gnomadg >> ${STEP5_DIR}/${file_label}_variant_annotation.log 2>&1; \
			echo "chrom\tpos\tref\talt\tgenotype\tgenotype_qual\tread_depth\tallele_depth\tvariant_allele_frac\tgenotype_likelihood\tVEP_Allele\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature_type\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tDISTANCE\tSTRAND\tFLAGS\tSYMBOL_SOURCE\tHGNC_ID\tCANONICAL\tMANE_SELECT\tMANE_PLUS_CLINICAL\tCCDS\tENSP\tRefSeq\tSIFT\tPolyPhen\tHGVS_OFFSET\tAF\tAFR_AF\tAMR_AF\tEAS_AF\tEUR_AF\tSAS_AF\tgnomADg_AF\tgnomADg_AFR_AF\tgnomADg_AMI_AF\tgnomADg_AMR_AF\tgnomADg_ASJ_AF\tgnomADg_EAS_AF\tgnomADg_FIN_AF\tgnomADg_MID_AF\tgnomADg_NFE_AF\tgnomADg_OTH_AF\tgnomADg_SAS_AF\tMAX_AF\tMAX_AF_POPS\tCLIN_SIG\tSOMATIC\tPHENO\tPUBMED" > ${STEP5_DIR}/${file_label}_PASS_NORM_PHASED_ANNOTATED_variants.tsv; \
			bcftools +split-vep ${STEP5_DIR}/${file_label}_PASS_NORM_PHASED_ANNOTATED_variants.vcf -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%GQ\t%DP\t%AD\t%VAF\t%PL]\t%CSQ\n' -d -A tab >> ${STEP5_DIR}/${file_label}_PASS_NORM_PHASED_ANNOTATED_variants.tsv; \
		fi; \
	fi
	@printf "\nFinished annotating small variants.\n";

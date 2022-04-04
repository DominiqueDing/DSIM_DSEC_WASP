#!/bin/bash

###################################################
### call Variants  ###############
###################################################
 ### adding GATK to PATH
export PATH=$PATH:/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/gatk-4.1.4.1/
# ### prepare reference fasta file
# cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference
# samtools faidx dsim-all-chromosome-r2.01.fasta
# java -jar /scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/picard/build/libs/picard.jar CreateSequenceDictionary \
# R=dsim-all-chromosome-r2.01.fasta \
# O=dsim-all-chromosome-r2.01.dict

cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/call_variants

 ### call variants round 1
echo "START CALLING VARIANTS ROUND 1..."
gatk --java-options "-Xmx4g" HaplotypeCaller \
-R /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference/dsim-all-chromosome-r2.01.fasta \
-I /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/genome_sequencing/dsim49_map_to_REFdsim-r2.01_sorted.dedup.sorted2.bam \
-O dsim49_1.vcf
### seperate snp and indel in vcf
gatk SelectVariants \
-V dsim49_1.vcf \
-select-type SNP \
-O dsim49_1_snp.vcf

gatk SelectVariants \
-V dsim49_1.vcf \
-select-type INDEL \
-O dsim49_1_indel.vcf

### hard-filtering
# for snp
gatk VariantFiltration \
-V dsim49_1_snp.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O dsim49_1_snp_hard-marked.vcf
# for indel
gatk VariantFiltration \
-V dsim49_1_indel.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "FS > 200.0" --filter-name "FS200" \
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
-O dsim49_1_indel_hard-marked.vcf

### select hard-filtered variants to serve as "true set" in the following process
# for snp
gatk SelectVariants \
-R /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference/dsim-all-chromosome-r2.01.fasta \
-V dsim49_1_snp_hard-marked.vcf \
-O dsim49_1_snp_hard-filtered.vcf \
--exclude-filtered
# for indel
gatk SelectVariants \
-R /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference/dsim-all-chromosome-r2.01.fasta \
-V dsim49_1_indel_hard-marked.vcf \
-O dsim49_1_indel_hard-filtered.vcf \
--exclude-filtered

### BSQR
gatk BaseRecalibrator \
-R /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference/dsim-all-chromosome-r2.01.fasta \
-I /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/genome_sequencing/dsim49_map_to_REFdsim-r2.01_sorted.dedup.sorted2.bam \
-known-sites dsim49_1_snp_hard-marked.vcf \
-O bsqr_1_recal_data.table
gatk ApplyBQSR \
-R /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference/dsim-all-chromosome-r2.01.fasta \
-I /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/genome_sequencing/dsim49_map_to_REFdsim-r2.01_sorted.dedup.sorted2.bam \
--bqsr-recal-file bsqr_1_recal_data.table \
-O bsqr_1_recal-ed.bam

##########################
### loops ################
##########################

declare -i count=2
while [[ count -le $1 ]]; do
	echo "STAR CALLING VARIANT ROUND ${count}..."
	### call variants with Haplotypecaller
	gatk --java-options "-Xmx4g" HaplotypeCaller \
-R /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference/dsim-all-chromosome-r2.01.fasta \
-I bsqr_$((count-1))_recal-ed.bam \
-O dsim49_"${count}".vcf
	### seperate snp and indel in vcf
	gatk SelectVariants \
-V dsim49_"${count}".vcf \
-select-type SNP \
-O dsim49_"${count}"_snp.vcf
	gatk SelectVariants \
-V dsim49_"${count}".vcf \
-select-type INDEL \
-O dsim49_"${count}"_indel.vcf
	### hard-filtering for snp
	gatk VariantFiltration \
-V dsim49_"${count}"_snp.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O dsim49_"${count}"_snp_hard-marked.vcf
	### hard-filtering for indel
	gatk VariantFiltration \
-V dsim49_"${count}"_indel.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "FS > 200.0" --filter-name "FS200" \
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
-O dsim49_"${count}"_indel_hard-marked.vcf
	### select hard-filtered variants to serve as "true set" in the following process
# for snp
	gatk SelectVariants \
-R /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference/dsim-all-chromosome-r2.01.fasta \
-V dsim49_"${count}"_snp_hard-marked.vcf \
-O dsim49_"${count}"_snp_hard-filtered.vcf \
--exclude-filtered;
# for indel
	gatk SelectVariants \
-R /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference/dsim-all-chromosome-r2.01.fasta \
-V dsim49_"${count}"_indel_hard-marked.vcf \
-O dsim49_"${count}"_indel_hard-filtered.vcf \
--exclude-filtered
	### BSQR
	gatk BaseRecalibrator \
-R /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference/dsim-all-chromosome-r2.01.fasta \
-I bsqr_$((count-1))_recal-ed.bam \
-known-sites dsim49_"${count}"_snp_hard-marked.vcf \
-O bsqr_"${count}"_recal_data.table
	gatk ApplyBQSR \
-R /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference/dsim-all-chromosome-r2.01.fasta \
-I bsqr_$((count-1))_recal-ed.bam \
--bqsr-recal-file bsqr_"${count}"_recal_data.table \
-O bsqr_"${count}"_recal-ed.bam
	((count++))
done






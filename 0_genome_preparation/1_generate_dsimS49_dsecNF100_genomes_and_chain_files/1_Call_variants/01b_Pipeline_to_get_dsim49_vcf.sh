#################
## get_assemly ##
#################

# convert sra file to fastq file:
fastq-dump SRR3585445.1

# quality check sequencing reads:
fastqc SRR3585445.1.fastq

###################################################
### map reads to reference genome:  ###############
###################################################

# using BWA so that read group information can be added during the alignment stage, without requiring a separate step with Picard

# create index with BWA
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference
bwa index -p dsim-r2.01 dsim-all-chromosome-r2.01.fasta
# mapping to dsim reference genome (bwa default setting)
bwa mem -M -t 1 -R '@RG\tID:C4K51ACXX.5\tPU:C4K51ACXX.5.sz49\tSM:sz49\tPL:ILLUMINA\tLB:sz49' dsim-r2.01 /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/genome_sequencing/SRR3585445.1.fastq > /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/genome_sequencing/dsim49_map_to_REFdsim-r2.01.sam

###################################################
### sorting and indexing aligned file  ###############
###################################################
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/genome_sequencing/
java -jar /scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/picard/build/libs/picard.jar SortSam \
I=dsim49_map_to_REFdsim-r2.01.sam \
O=dsim49_map_to_REFdsim-r2.01_sorted.bam \
SORT_ORDER=coordinate \
CREATE_INDEX=true

###################################################
### Alignment metrics  ###############
###################################################
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/genome_sequencing/
java -jar /scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/picard/build/libs/picard.jar CollectAlignmentSummaryMetrics \
I=dsim49_map_to_REFdsim-r2.01_sorted.bam \
R=/scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference/dsim-all-chromosome-r2.01.fasta \
METRIC_ACCUMULATION_LEVEL=SAMPLE \
METRIC_ACCUMULATION_LEVEL=READ_GROUP \
O=dsim49_map_to_REFdsim-r2.01.alignment_metrics.txt

###################################################
### mark duplicates  ###############
###################################################
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/genome_sequencing/
java -jar /scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/picard/build/libs/picard.jar MarkDuplicates \
TMP_DIR=tmp \
I=dsim49_map_to_REFdsim-r2.01_sorted.bam \
O=dsim49_map_to_REFdsim-r2.01_sorted.dedup.bam \
METRICS_FILE=dsim49_map_to_REFdsim-r2.01_sorted.dedup_metrics.txt \
REMOVE_DUPLICATES=false \
TAGGING_POLICY=All
# sort and index bam file again:
java -jar /scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/picard/build/libs/picard.jar SortSam \
I=dsim49_map_to_REFdsim-r2.01_sorted.dedup.bam \
O=dsim49_map_to_REFdsim-r2.01_sorted.dedup.sorted2.bam \
SORT_ORDER=coordinate \
CREATE_INDEX=true

###################################################
### call Variants  ###############
###################################################
### adding GATK to PATH
export PATH=$PATH:/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/gatk-4.1.4.1/
### prepare reference fasta file
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference
samtools faidx dsim-all-chromosome-r2.01.fasta
java -jar /scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/picard/build/libs/picard.jar CreateSequenceDictionary \
R=dsim-all-chromosome-r2.01.fasta \
O=dsim-all-chromosome-r2.01.dict

### calling variant round1
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/call_variants

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
-O dsim49_1_snp_hard-filtered.vcf
# for indel
gatk VariantFiltration \
-V dsim49_1_indel.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "FS > 200.0" --filter-name "FS200" \
-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
-O dsim49_1_indel_hard-filtered.vcf

### select hard-filtered variants to serve as "true set" in the following process
# for snp
gatk SelectVariants \
-R /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference/dsim-all-chromosome-r2.01.fasta \
-V dsim49_1_snp_hard-filtered.vcf \
-O dsim49_1_snp_custom-true.vcf \
--exclude-filtered
# for indel
gatk SelectVariants \
-R /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference/dsim-all-chromosome-r2.01.fasta \
-V dsim49_1_indel_hard-filtered.vcf \
-O dsim49_1_indel_custom-true.vcf \
--exclude-filtered

### Calculate VQSLOD tranches
# for snp
gatk --java-options "-Xmx3g -Xms3g" VariantRecalibrator \
-V dsim49_1.vcf \
--trust-all-polymorphic \
-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
-mode SNP \
--max-gaussians 6 \
--resource:custom,known=false,training=true,truth=true,prior=10.0 dsim49_1_snp_custom-true.vcf \
-O dsim49_1_snp.recal \
--tranches-file dsim49_1_snp.tranches
#-resource hapmap,known=false,training=true,truth=true,prior=15:hapmap_3.3.hg38.vcf.gz \
#-resource omni,known=false,training=true,truth=true,prior=12:1000G_omni2.5.hg38.vcf.gz \
#-resource 1000G,known=false,training=true,truth=false,prior=10:1000G_phase1.snps.high_confidence.hg38.vcf.gz \
#-resource dbsnp,known=true,training=false,truth=false,prior=7:Homo_sapiens_assembly38.dbsnp138.vcf \
# for indels
gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
-V dsim49_1.vcf \
--trust-all-polymorphic \
-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
-mode INDEL \
--max-gaussians 4 \
--resource:custom,known=false,training=true,truth=true,prior=10.0 dsim49_1_indel_custom-true.vcf \
-O dsim49_1_indel.recal \
--tranches-file dsim49_1_indel.tranches
#-resource mills,known=false,training=true,truth=true,prior=12:Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
#-resource axiomPoly,known=false,training=true,truth=false,prior=10:Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
#-resource dbsnp,known=true,training=false,truth=false,prior=2:Homo_sapiens_assembly38.dbsnp138.vcf \

###  Filter on VQSLOD using ApplyVQSR
# for indel
gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
-V dsim49_1.vcf \
--recal-file dsim49_1_indel.recal \
--tranches-file dsim49_1_indel.tranches \
--truth-sensitivity-filter-level 99.7 \
--create-output-variant-index true \
-mode INDEL \
-O dsim49_1_indel-recalibrated.vcf
# for snp
gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
-V dsim49_1_indel-recalibrated.vcf \
--recal-file dsim49_1_snp.recal \
--tranches-file dsim49_1_snp.tranches \
--truth-sensitivity-filter-level 99.7 \
--create-output-variant-index true \
-mode SNP \
-O dsim49_1_indel-snp-recalibrated.vcf

### select passed variants
gatk SelectVariants \
-R /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference/dsim-all-chromosome-r2.01.fasta \
-V dsim49_1_indel-snp-recalibrated.vcf \
-O dsim49_2.vcf \
--exclude-filtered






















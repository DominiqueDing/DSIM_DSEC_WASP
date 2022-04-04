#!/bin/bash

# mapping sequence: 1-->4-->2-->3

#  # indexing v2 reference genomes with bowtie: 

# cd /scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map/dsecNF100v2_genome/
# cp /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff_2/liftover/dsecNF100v2_genome/dsecNF100v2.* ./
# bowtie-build dsecNF100v2.fasta sec_bowtie

# cd /scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map/dsim49v2_genome/
# cp /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim49v2_genome/dsim49v2.* ./
# bowtie-build dsim49v2.fasta sim_bowtie

# cd /scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map/dsimREF_genome/
# bowtie-build dsim-all-chromosome-r2.01.fasta simREF_bowtie

# cd /scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map/dsecREF_genome/
# bowtie-build dsec-all-chromosome-r1.3.fasta secREF_bowtie


BAR=A
LANE='s_1'


PREFIX_3='maptodsimREFr2'
PREFIX_2='maptodsim49v2r2'
PREFIX_4='maptodsecREFr1'
PREFIX_1='maptodsecNF100v2r1'


OP_DIR_3='/scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map/dsimREF_genome'
OP_DIR_2='/scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map/dsim49v2_genome'
OP_DIR_4='/scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map/dsecREF_genome'
OP_DIR_1='/scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map/dsecNF100v2_genome'

OUT_DIR='/scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map/map_map3_v2_genomes'

UN_3='simREF_map2_results/unaligned_reads'
AL_3='simREF_map2_results/aligned_reads'
UN_2='sim_map2_results/unaligned_reads'
AL_2='sim_map2_results/aligned_reads'
UN_4='secREF_map1_results/unaligned_reads'
AL_4='secREF_map1_results/aligned_reads'
UN_1='sec_map1_results/unaligned_reads'
AL_1='sec_map1_results/aligned_reads'

SAM_3='simREF_map2_results/sam_files'
SAM_2='sim_map2_results/sam_files'
SAM_4='secREF_map1_results/sam_files'
SAM_1='sec_map1_results/sam_files'

BOWTIE_INDEX_3='simREF_bowtie'
BOWTIE_INDEX_2='sim_bowtie'
BOWTIE_INDEX_4='secREF_bowtie'
BOWTIE_INDEX_1='sec_bowtie'

IN_DIR_FQ='/scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map'


########################################
### map to secNF100.2v2 round1 ###########
########################################
cd ${OP_DIR_1}

for i in {1..12}; do \
bowtie \
-v 0 -k 3 --best --strata -m 3 \
--un "${OUT_DIR}"/"${UN_1}"/"${PREFIX_1}".un1."${BAR}${i}"."${LANE}".fq \
--max "${OUT_DIR}"/"${UN_1}"/"${PREFIX_1}".max1."${BAR}${i}"."${LANE}".fq \
--al "${OUT_DIR}"/"${AL_1}"/"${PREFIX_1}".al1."${BAR}${i}"."${LANE}".fq \
--sam --no-unal \
"${BOWTIE_INDEX_1}" \
"${IN_DIR_FQ}"/"${LANE}"/SLX-17988."${BAR}${i}".HL7FNDRXX."${LANE}".r_1_trimmed.fq \
"${OUT_DIR}"/"${SAM_1}"/"${PREFIX_1}"."${BAR}${i}"."${LANE}".sam; done

######################################
### map to secREF round1 #############
######################################
cd ${OP_DIR_4}

for i in {1..12}; do \
bowtie \
-v 0 -k 3 --best --strata -m 3 \
--un "${OUT_DIR}"/"${UN_4}"/"${PREFIX_4}".un1."${BAR}${i}"."${LANE}".fq \
--max "${OUT_DIR}"/"${UN_4}"/"${PREFIX_4}".max1."${BAR}${i}"."${LANE}".fq \
--al "${OUT_DIR}"/"${AL_4}"/"${PREFIX_4}".al1."${BAR}${i}"."${LANE}".fq \
--sam --no-unal \
"${BOWTIE_INDEX_4}" \
"${OUT_DIR}"/"${UN_1}"/"${PREFIX_1}".un1."${BAR}${i}"."${LANE}".fq \
"${OUT_DIR}"/"${SAM_4}"/"${PREFIX_4}"."${BAR}${i}"."${LANE}".sam; done






#####################################
### map to sim49v2 round2 #############
#####################################
cd ${OP_DIR_2}

for i in {1..12}; do \
bowtie \
-v 0 -k 3 --best --strata -m 3 \
--un "${OUT_DIR}"/"${UN_2}"/"${PREFIX_2}".un2."${BAR}${i}"."${LANE}".fq \
--max "${OUT_DIR}"/"${UN_2}"/"${PREFIX_2}".max2."${BAR}${i}"."${LANE}".fq \
--al "${OUT_DIR}"/"${AL_2}"/"${PREFIX_2}".al2."${BAR}${i}"."${LANE}".fq \
--sam --no-unal \
"${BOWTIE_INDEX_2}" \
"${OUT_DIR}"/"${UN_4}"/"${PREFIX_4}".un1."${BAR}${i}"."${LANE}".fq \
"${OUT_DIR}"/"${SAM_2}"/"${PREFIX_2}"."${BAR}${i}"."${LANE}".sam; done

#####################################
### map to simREF round2 ############
#####################################
cd ${OP_DIR_3}

for i in {1..12}; do \
bowtie \
-v 0 -k 3 --best --strata -m 3 \
--un "${OUT_DIR}"/"${UN_3}"/"${PREFIX_3}".un2."${BAR}${i}"."${LANE}".fq \
--max "${OUT_DIR}"/"${UN_3}"/"${PREFIX_3}".max2."${BAR}${i}"."${LANE}".fq \
--al "${OUT_DIR}"/"${AL_3}"/"${PREFIX_3}".al2."${BAR}${i}"."${LANE}".fq \
--sam --no-unal \
"${BOWTIE_INDEX_3}" \
"${OUT_DIR}"/"${UN_2}"/"${PREFIX_2}".un2."${BAR}${i}"."${LANE}".fq \
"${OUT_DIR}"/"${SAM_3}"/"${PREFIX_3}"."${BAR}${i}"."${LANE}".sam; done


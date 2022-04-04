#!/bin/bash

# extract uniquely-mapped reads:
MASTER_DIR="/scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map/map_dsimREF_dsecREF_salvage"
OUT_DIR="/scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map/map_dsimREF_dsecREF_salvage/Annotating/MAP3_annot"

SEC_DIR="sec_map2_results"
SECREF_DIR="secREF_map2_results"
SIM_DIR="sim_map2_results"
SIMREF_DIR="simREF_map2_results"

SEC_PRE="maptodsecNF100r2"
SECREF_PRE="maptodsecREFr2"
SIM_PRE="maptodsim49r2"
SIMREF_PRE="maptodsimREFr2"



for BAR in A B C D; do
	for i in {1..12}; do
		for lane in s_1 s_2; do
			grep "XM:i:2" "${MASTER_DIR}"/"${SEC_DIR}"/sam_files/"${SEC_PRE}"."${BAR}${i}"."${lane}".sam > "${OUT_DIR}"/"${SEC_PRE}"."${BAR}${i}"."${lane}".uniq.sam
			sam2bed <  "${OUT_DIR}"/"${SEC_PRE}"."${BAR}${i}"."${lane}".uniq.sam > "${OUT_DIR}"/"${SEC_PRE}"."${BAR}${i}"."${lane}".uniq.bed
		done
	done
done


for BAR in A B C D; do
	for i in {1..12}; do
		for lane in s_1 s_2; do
			grep "XM:i:2" "${MASTER_DIR}"/"${SECREF_DIR}"/sam_files/"${SECREF_PRE}"."${BAR}${i}"."${lane}".sam > "${OUT_DIR}"/"${SECREF_PRE}"."${BAR}${i}"."${lane}".uniq.sam
			sam2bed <  "${OUT_DIR}"/"${SECREF_PRE}"."${BAR}${i}"."${lane}".uniq.sam > "${OUT_DIR}"/"${SECREF_PRE}"."${BAR}${i}"."${lane}".uniq.bed
		done
	done
done


for BAR in A B C D; do
	for i in {1..12}; do
		for lane in s_1 s_2; do
			grep "XM:i:2" "${MASTER_DIR}"/"${SIM_DIR}"/sam_files/"${SIM_PRE}"."${BAR}${i}"."${lane}".sam > "${OUT_DIR}"/"${SIM_PRE}"."${BAR}${i}"."${lane}".uniq.sam
			sam2bed <  "${OUT_DIR}"/"${SIM_PRE}"."${BAR}${i}"."${lane}".uniq.sam > "${OUT_DIR}"/"${SIM_PRE}"."${BAR}${i}"."${lane}".uniq.bed
		done
	done
done


for BAR in A B C D; do
	for i in {1..12}; do
		for lane in s_1 s_2; do
			grep "XM:i:2" "${MASTER_DIR}"/"${SIMREF_DIR}"/sam_files/"${SIMREF_PRE}"."${BAR}${i}"."${lane}".sam > "${OUT_DIR}"/"${SIMREF_PRE}"."${BAR}${i}"."${lane}".uniq.sam
			sam2bed <  "${OUT_DIR}"/"${SIMREF_PRE}"."${BAR}${i}"."${lane}".uniq.sam > "${OUT_DIR}"/"${SIMREF_PRE}"."${BAR}${i}"."${lane}".uniq.bed
		done
	done
done
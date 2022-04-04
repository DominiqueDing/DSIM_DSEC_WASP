#!/bin/bash


############################################
########## trim with trim-galore: ##########
############################################
cd /scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map
source ~/use_python-3.8.2.sh
for i in `cat lib.name.2.txt`; do \
/home/sd768/download_prog/TrimGalore-0.6.0/trim_galore --cores 4 --quality 30 \
--phred33 --fastqc --illumina --length 50 --output_dir s_1/ \
--path_to_cutadapt /home/sd768/.local/bin/cutadapt /home/public/Dsim_Dsec_Wasp_raw_sequencing/SLX-17988."${i}".HL7FNDRXX.s_1.r_1.fq.gz; \
/home/sd768/download_prog/TrimGalore-0.6.0/trim_galore --cores 4 --quality 30 \
--phred33 --fastqc --illumina --length 50 --output_dir s_2/ \
--path_to_cutadapt /home/sd768/.local/bin/cutadapt /home/public/Dsim_Dsec_Wasp_raw_sequencing/SLX-17988."${i}".HL7FNDRXX.s_2.r_1.fq.gz; done
# catenate lanes for eachh lib:
# mkdir cat_s1s2_trimmed 
# for i in {1..12}; do \
# cat s_1/SLX-17988.E${i}.HL7FNDRXX.s_1.r_1_trimmed.fq.gz s_2/SLX-17988.E${i}.HL7FNDRXX.s_2.r_1_trimmed.fq.gz > \
# cat_s1s2_trimmed/E${i}_trimmed.fq.gz; done


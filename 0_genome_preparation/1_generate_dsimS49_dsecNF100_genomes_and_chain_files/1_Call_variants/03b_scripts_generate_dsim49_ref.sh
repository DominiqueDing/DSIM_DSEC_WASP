######################################
## generate alternative ref dsim49 ###
######################################
### adding GATK to PATH
export PATH=$PATH:/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/gatk-4.1.4.1/
###
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim49_genome

### generate new dsim49 reference fasta:
gatk FastaAlternateReferenceMaker \
-R /scratch/Dominique/dsim_dsec_wasp_proj/dsim_reference/dsim-all-chromosome-r2.01.fasta \
-O dsim49.fasta \
-V dsim49_24.vcf




######################################
## generate alternative ref dsecNF100 ###
######################################
### adding GATK to PATH
export PATH=$PATH:/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/gatk-4.1.4.1/
###
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff_2/liftover/dsecNF100_genome

### generate new dsecNF100 reference fasta:
gatk FastaAlternateReferenceMaker \
-R /scratch/Dominique/dsim_dsec_wasp_proj/dsec_reference/dsec-all-chromosome-r1.3.fasta \
-O dsecNF100.fasta \
-V dsecNF100_20.vcf



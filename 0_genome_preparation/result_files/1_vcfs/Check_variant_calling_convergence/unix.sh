# get GQ values from all vcf files:
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff_2/call_variants
for i in {1..20}; do more +/#CHROM dsecNF100_${i}.vcf | awk '{print $NF}' | awk 'BEGIN{FS=":"}{print $4}' | sed -e 's/$/\t'${i}'/' > dsecNF100.round.${i}.BQSR.GQ.txt;done

cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/call_variants
for i in {1..24}; do more +/#CHROM dsim49_${i}.vcf | awk '{print $NF}' | awk 'BEGIN{FS=":"}{print $4}' | sed -e 's/$/\t'${i}'/' > dsim49.round.${i}.BQSR.GQ.txt;done
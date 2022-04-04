
# go to operating directory
cd /Volumes/Treasure_Island/PROJECT_wasp_dsim_dsec/Liftover/R_get_constitutive
# get mRNA-gene from gff
grep "\tmRNA\t" ../011dsec_r1.3/dsec-all-no-analysis-r1.3.gff| awk '{print $NF}'|sed 's/ID=//'|sed -e 's/;.*Parent=/'$'\t''/'|awk 'BEGIN{FS=";"}{print $1}' > dsecREF.mRNA.gene.tab
# get exon-mRNA from gff
grep "\texon\t" ../011dsec_r1.3/dsec-all-no-analysis-r1.3.gff| sed 's/'$'\t''/:/g'|sed 's/Parent=/'$'\t''/' > dsecREF.exon.mRNA.tab 


# get mRNA-gene from gff
grep "\tmRNA\t" ../012_dsim_r2.01/dsim-all-r2.01.gff| awk '{print $NF}'|sed 's/ID=//'|sed -e 's/;.*Parent=/'$'\t''/'|awk 'BEGIN{FS=";"}{print $1}' > dsimREF.mRNA.gene.tab
# get exon-mRNA from gff
grep "\texon\t" ../012_dsim_r2.01/dsim-all-r2.01.gff| sed 's/'$'\t''/:/g'|awk 'BEGIN{FS=";"}{print $1}'|sed 's/Parent=/'$'\t''/'|sed 's/,/'$'\t''/g' > dsimREF.exon.mRNA.tab 
awk '{print NF}' dsimREF.exon.mRNA.tab|sort -nu|tail -n 1 # find out highest column counts: 45



# g
sed -i '' 's/:/'$'\t''/g' dsimREF.constitutive.exons.txt
sed -i '' 's/:/'$'\t''/g' dsecREF.constitutive.exons.txt
# grep done on multivac "/scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map/Get_constitutive_exons", too CPU demanding...
grep -f dsimREF.constitutive.exons.txt dsim-all-r2.01.gtf > dsimREF.constitutive.exons.gtf
grep -f dsecREF.constitutive.exons.txt dsec-all-r1.3.gtf > dsecREF.constitutive.exons.gtf

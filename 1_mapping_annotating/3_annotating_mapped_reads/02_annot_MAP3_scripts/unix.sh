# generate scripts for parallel annotating
for bar in A B C D; \
do for tag in dsecNF100 dsecREF dsim49 dsimREF; \
do sed 's/annot_eg/MAP3_annot/g' for.mapto${tag}r2.sh | sed 's/BAR="."/BAR="'${bar}'"/g' |sed 's/for i in .;/for i in {1..12};/g' > annot_MAP3_scripts/${bar}.for.mapto${tag}r2.sh;done;done 

# gather counts from *OnSecREF.OnSimREF.const_counts files:
for gen in dsecNF100 dsecREF dsim49 dsimREF; \
do rm ${gen}.conserved.constitutive.counts; for bar in A B C D; \
do for i in {1..12}; \
do for lane in s_1 s_2; \
do sed 's/$/'"\t"${gen}"\t"${bar}${i}'/g' mapto${gen}/mapto${gen}r2.${bar}${i}.${lane}.uniq.OnSecREF.OnSimREF.const_counts >> conserved.constitutive.counts; \
done; done; done; done

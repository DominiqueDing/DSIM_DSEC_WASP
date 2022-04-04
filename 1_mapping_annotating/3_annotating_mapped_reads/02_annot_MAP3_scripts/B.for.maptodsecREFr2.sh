chain_dir="/scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map/map_dsimREF_dsecREF_salvage/Annotating/chain_files"
genome_dir="/scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map/map_dsimREF_dsecREF_salvage/Annotating/genome_files"
opt_dir_home="/scratch/Dominique/dsim_dsec_wasp_proj/trim_amd_map/map_dsimREF_dsecREF_salvage/Annotating/MAP3_annot"
opt_dir="maptodsecREF"

map_prefix="maptodsecREFr2"
BAR="B"
lane="s_1"




for i in {1..12}; do \
#########################
## map to dsecREF.const #
#########################
rm "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.bed;\
htseq-count -f sam --stranded=no -o "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREF.htseq-out.sam \
"${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.sam \
"${genome_dir}"/dsecREF.constitutive.exons.gtf > "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREF.const_counts;\
#rm "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.sam;\
grep "XF:Z:FBgn" "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREF.htseq-out.sam > "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.sam;\
grep "XF:Z:__ambiguous" "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREF.htseq-out.sam > "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConstAmbi.sam;\
rm "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREF.htseq-out.sam;\


#########################
## liftover to dsimREF ##
#########################
sam2bed <  "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.sam > "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.bed;\
rm "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.sam;\
liftOver -bedPlus=3 "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.bed \
"${chain_dir}"/droSec1ToDroSim2.over.Modi.chain \
"${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREF.bed \
"${opt_dir_home}/${opt_dir}"/Unlift;\
rm "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.bed;\

bedToBam -i "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREF.bed \
-g "${genome_dir}"/dsimREF.chrom.sizes | samtools view > "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREF.incomplete.sam;\

./GetStringforSAMfromBED "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREF.bed \
"${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREF.incomplete.sam \
"${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREF.sam; \

# while IFS= read LINE; \
# do SEQ=`echo "${LINE}"|awk '{print $1}'`;\
# F10=`grep "${SEQ}" "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREF.bed|awk '{print $12}'`;\
# F11=`grep "${SEQ}" "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREF.bed|awk '{print $13}'`; \
# echo "${LINE}"|awk -F"\t" 'BEGIN{OFS="\t";}{$10="'${F10}'";}{$11="'${F11}'";}{print;}' >> "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREF.sam;\
# done < "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREF.incomplete.sam;\

rm "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREF.incomplete.sam;\
rm "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREF.bed;\


#########################
## map to dsimREF.const #
#########################
htseq-count -f sam --stranded=no -o "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREFConst.htseq-out.sam \
"${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREF.sam \
"${genome_dir}"/dsimREF.constitutive.exons.gtf > "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREF.OnSimREF.const_counts;\
rm "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREF.sam;\
grep "XF:Z:FBgn" "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREFConst.htseq-out.sam > "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREFConst.sam;\
grep "XF:Z:__ambiguous" "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREFConst.htseq-out.sam > "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConstOnSimREFConstAmbi.sam;\
rm "${opt_dir_home}/${opt_dir}"/"${map_prefix}"."${BAR}${i}"."${lane}".uniq.OnSecREFConst.OnSimREFConst.htseq-out.sam;\

done

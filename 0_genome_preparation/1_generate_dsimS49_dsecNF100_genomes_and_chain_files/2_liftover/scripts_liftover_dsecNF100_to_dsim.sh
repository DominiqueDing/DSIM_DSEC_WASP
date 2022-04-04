# make directories for dsim and dsecNF100 assembli:
mkdir dsim
mkdir dsecNF100
cp path_to_dsim_ref/dsim-all-chromosome-r2.02.fasta dsim/
cp /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff/liftover/dsecNF100_1_genome/dsecNF100_1_ref.fasta dsecNF100/
cat dsecNF100_1_ref.fasta | sed 's/>.* />dsecNF100_/' | sed 's/:/ len:/' > my_dsecNF100_1_ref.fasta
# create twobit file for dsim_ref:
cd dsim
/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/ucsc_utilities/faToTwoBit dsim-all-chromosome-r2.02.fasta dsim.2bit
/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/ucsc_utilities/twoBitInfo dsim.2bit stdout | sort -k2rn > dsim.chrom.sizes
# create twobit file for dsecNF100:
cd dsecNF100
/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/ucsc_utilities/faToTwoBit my_dsecNF100_1_ref.fasta my_dsecNF100.2bit
/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/ucsc_utilities/twoBitInfo my_dsecNF100.2bit stdout | sort -k2rn > my_dsecNF100.chrom.sizes
# split dsecNF100:
mkdir split_fas
cd split_fas
faSplit gap ../my_dsecNF100_1_ref.fasta 50000 dsecNF100- -lift=dsecNF100-split.lft
# blat against dsim_ref:
mkdir pls_0
cd split_fas
for i in *.fa; do blat /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff/liftover/dsim/dsim.2bit ${i} ../pls_0/"OLD."`echo ${i} | sed 's/.fa/.psl/'`; done
# liftUp change the coordinate system:
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff/liftover/dsecNF100/
mkdir pls_1
cd pls_0
for i in *.psl; do liftUp -pslQ ../pls_1/`echo ${i} | sed 's/OLD.//'` /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff/liftover/dsecNF100/split_fas/dsecNF100-split.lft warn ${i}; done
# translate pls to chain file:
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff/liftover/dsecNF100
mkdir dsim_to_dsecNF100_chain
cd pls_1
for i in *.psl; do axtChain -linearGap=medium -psl ${i} \
/scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff/liftover/dsim/dsim.2bit \
/scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff/liftover/dsecNF100/my_dsecNF100.2bit ../dsim_to_dsecNF100_chain/`echo ${i} | sed 's/psl/chain/'`; done
# Combine and sort chain files:
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff/liftover/dsecNF100
chainMergeSort dsim_to_dsecNF100_chain/*.chain | chainSplit dsim_to_dsecNF100_merged_chain stdin
# Concat and sort the chains:
cat dsim_to_dsecNF100_merged_chain/*.chain > dsim_to_dsecNF100_merged_chain/all.chain
chainSort dsim_to_dsecNF100_merged_chain/all.chain dsim_to_dsecNF100_merged_chain/all.sorted.chain
# Netting: identify alignable regions from chains:
mkdir dsim_to_dsecNF100_net
chainNet dsim_to_dsecNF100_merged_chain/all.sorted.chain \
/scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff/liftover/dsim/dsim.chrom.sizes \
/scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff/liftover/dsecNF100/my_dsecNF100.chrom.sizes \
dsim_to_dsecNF100_net/all.net /dev/null
#  Finally, select the right alignable regions using the nets, creating a "liftOver" file:
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff/liftover/dsecNF100
mkdir DsimTodsecNF100_liftover_chain_file
netChainSubset dsim_to_dsecNF100_net/all.net dsim_to_dsecNF100_merged_chain/all.chain DsimTodsecNF100_liftover_chain_file/DsimTodsecNF100_liftover.chain




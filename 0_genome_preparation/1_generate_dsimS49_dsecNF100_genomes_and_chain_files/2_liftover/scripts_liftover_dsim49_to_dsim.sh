# make directories for dsim and dsim49 assembli:
mkdir dsim
mkdir dsim49
cp path_to_dsim_ref/dsim-all-chromosome-r2.02.fasta dsim/
cp /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim49_1_genome/dsim49_1_ref.fasta dsim49/
cat dsim49_1_ref.fasta | sed 's/>.* />dsim49_/' | sed 's/:/ len:/' > my_dsim49_1_ref.fasta
# create twobit file for dsim_ref:
cd dsim
/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/ucsc_utilities/faToTwoBit dsim-all-chromosome-r2.02.fasta dsim.2bit
/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/ucsc_utilities/twoBitInfo dsim.2bit stdout | sort -k2rn > dsim.chrom.sizes
# create twobit file for dsim49:
cd dsim49
/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/ucsc_utilities/faToTwoBit my_dsim49_1_ref.fasta my_dsim49.2bit
/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/ucsc_utilities/twoBitInfo my_dsim49.2bit stdout | sort -k2rn > my_dsim49.chrom.sizes
# split dsim49:
mkdir split_fas
cd split_fas
faSplit gap ../my_dsim49_1_ref.fasta 50000 dsim49- -lift=dsim49-split.lft
# blat against dsim_ref:
mkdir pls_0
cd split_fas
for i in *.fa; do blat /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim/dsim.2bit ${i} ../pls_0/"OLD."`echo ${i} | sed 's/.fa/.psl/'`; done
# liftUp change the coordinate system:
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim49/
mkdir pls_1
cd pls_0
for i in *.psl; do liftUp -pslQ ../pls_1/`echo ${i} | sed 's/OLD.//'` /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim49/split_fas/dsim49-split.lft warn ${i}; done
# translate pls to chain file:
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim49
mkdir dsim_to_dsim49_chain
cd pls_1
for i in *.psl; do axtChain -linearGap=medium -psl ${i} \
/scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim/dsim.2bit \
/scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim49/my_dsim49.2bit ../dsim_to_dsim49_chain/`echo ${i} | sed 's/psl/chain/'`; done
# Combine and sort chain files:
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim49
chainMergeSort dsim_to_dsim49_chain/*.chain | chainSplit dsim_to_dsim49_merged_chain stdin
# Concat and sort the chains:
cat dsim_to_dsim49_merged_chain/*.chain > dsim_to_dsim49_merged_chain/all.chain
chainSort dsim_to_dsim49_merged_chain/all.chain dsim_to_dsim49_merged_chain/all.sorted.chain
# Netting: identify alignable regions from chains:
mkdir dsim_to_dsim49_net
chainNet dsim_to_dsim49_merged_chain/all.sorted.chain \
/scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim/dsim.chrom.sizes \
/scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim49/my_dsim49.chrom.sizes \
dsim_to_dsim49_net/all.net /dev/null
#  Finally, select the right alignable regions using the nets, creating a "liftOver" file:
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim49
mkdir DsimToDsim49_liftover_chain_file
netChainSubset dsim_to_dsim49_net/all.net dsim_to_dsim49_merged_chain/all.chain DsimToDsim49_liftover_chain_file/DsimToDsim49_liftover.chain




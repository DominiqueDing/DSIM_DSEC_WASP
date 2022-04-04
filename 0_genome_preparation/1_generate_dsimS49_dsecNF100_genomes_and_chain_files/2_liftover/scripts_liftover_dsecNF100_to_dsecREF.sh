# make directories for dsim and dsecNF100 assembli:
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff_2/liftover
# create twobit file for dsec1.3:
cd dsec1.3_genome
/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/ucsc_utilities/faToTwoBit dsec-all-chromosome-r1.3.fasta dsec1.3.2bit
/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/ucsc_utilities/twoBitInfo dsec1.3.2bit stdout | sort -k2rn > dsec1.3.chrom.sizes
# create twobit file for dsecNF100:
cd ../dsecNF100_genome/
/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/ucsc_utilities/faToTwoBit dsecNF100.fasta dsecNF100.2bit
/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/ucsc_utilities/twoBitInfo dsecNF100.2bit stdout | sort -k2rn > dsecNF100.chrom.sizes

###########################################
## dsecNF100(OLD).to.dsec1.3(NEW).chain ###
###########################################
### Split NEW assembly
# split dsec1.3(NEW):
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff_2/liftover/dsec1.3_genome
mkdir split_fas
cd split_fas
faSplit size ../dsec-all-chromosome-r1.3.fasta 3000 dsecREF -lift=dsec1.3-split.lft -oneFile
# blat dsec1.3(NEW) against dsecNF100(Old):
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff_2/liftover/dsec1.3_genome
mkdir pls
cd split_fas
blat /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff_2/liftover/dsecNF100_genome/dsecNF100.2bit dsecREF.fa OLD.psl -tileSize=12 -minScore=100 -minIdentity=98 -fastMap
# liftUp change the coordinate system:
liftUp -pslQ dsecREF.psl dsec1.3-split.lft warn OLD.psl
# Chain together alignments 
axtChain -linearGap=medium -psl dsecREF.psl /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff_2/liftover/dsecNF100_genome/dsecNF100.2bit /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff_2/liftover/dsec1.3_genome/dsec1.3.2bit dsecREF.chain
# sort chain file
chainSort dsecREF.chain dsecREF.sorted.chain
# Netting: identify alignable regions from chains:
chainNet dsecREF.sorted.chain /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff_2/liftover/dsecNF100_genome/dsecNF100.chrom.sizes /scratch/Dominique/dsim_dsec_wasp_proj/dsecNF100_stuff_2/liftover/dsec1.3_genome/dsec1.3.chrom.sizes dsecREF.net /dev/null
# Finally, select the right alignable regions using the nets, creating a "liftOver" file:
netChainSubset dsecREF.net dsecREF.sorted.chain dsecNF100_To_dsecREF.over.chain


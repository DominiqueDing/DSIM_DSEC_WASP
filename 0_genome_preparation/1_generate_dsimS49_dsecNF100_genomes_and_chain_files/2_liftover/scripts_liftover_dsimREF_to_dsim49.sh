# make directories for dsim49(OLD) and dsimREF(NEW) assembli:
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover
# create twobit file for dsim2.02/dsimREF:
cd dsimREF_genome/
/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/ucsc_utilities/faToTwoBit dsim-all-chromosome-r2.02.fasta dsimREF.2bit
/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/ucsc_utilities/twoBitInfo dsimREF.2bit stdout | sort -k2rn > dsimREF.chrom.sizes
# create twobit file for dsim49(OLD):
cd ../dsim49_genome/
/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/ucsc_utilities/faToTwoBit dsim49.fasta dsim49.2bit
/scratch/Dominique/dsim_dsec_wasp_proj/downloaded_tools/ucsc_utilities/twoBitInfo dsim49.2bit stdout | sort -k2rn > dsim49.chrom.sizes

###########################################
## dsim2.02/dsimREF(OLD).to.dsim49(NEW).chain ###
###########################################
### Split NEW assembly
# split dsim49(NEW):
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim49_genome
mkdir split_fas
cd split_fas
faSplit size ../dsim49.fasta 3000 dsim49 -lift=dsim49-split.lft -oneFile
# blat dsim49(NEW) against dsimREF(OLD):
blat /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsimREF_genome/dsimREF.2bit dsim49.fa OLD.psl -tileSize=12 -minScore=100 -minIdentity=98 -fastMap
# liftUp change the coordinate system:
liftUp -pslQ dsim49.psl dsim49-split.lft warn OLD.psl
# Chain together alignments 
xtChain -linearGap=medium -psl dsim49.psl /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsimREF_genome/dsimREF.2bit /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim49_genome/dsim49.2bit dsim49.chain
# sort chain file
chainSort dsim49.chain dsim49.sorted.chain
# Netting: identify alignable regions from chains:
chainNet dsim49.sorted.chain /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsimREF_genome/dsimREF.chrom.sizes /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim49_genome/dsim49.chrom.sizes dsim49.net /dev/null
# Finally, select the right alignable regions using the nets, creating a "liftOver" file:
netChainSubset dsim49.net dsim49.sorted.chain dsimREF_To_dsim49.over.chain



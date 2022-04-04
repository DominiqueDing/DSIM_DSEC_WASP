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
## dsim49(OLD).to.dsim2.02/dsimREF(NEW).chain ###
###########################################
### Split NEW assembly
# split dsimREF(NEW):
cd /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsimREF_genome
mkdir split_fas
cd split_fas
faSplit size ../dsim-all-chromosome-r2.02.fasta 3000 dsimREF -lift=dsim2.02-split.lft -oneFile
# blat dsimREF(NEW) against dsim49(OLD):
blat /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim49_genome/dsim49.2bit dsimREF.fa OLD.psl -tileSize=12 -minScore=100 -minIdentity=98 -fastMap
# liftUp change the coordinate system:
liftUp -pslQ dsimREF.psl dsim2.02-split.lft warn OLD.psl
# Chain together alignments 
axtChain -linearGap=medium -psl dsimREF.psl /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim49_genome/dsim49.2bit /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsimREF_genome/dsimREF.2bit dsimREF.chain
# sort chain file
chainSort dsimREF.chain dsimREF.sorted.chain
# Netting: identify alignable regions from chains:
chainNet dsimREF.sorted.chain /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsim49_genome/dsim49.chrom.sizes /scratch/Dominique/dsim_dsec_wasp_proj/dsim49_stuff/liftover/dsimREF_genome/dsimREF.chrom.sizes dsimREF.net /dev/null
# Finally, select the right alignable regions using the nets, creating a "liftOver" file:
netChainSubset dsimREF.net dsimREF.sorted.chain dsim49_To_dsimREF.over.chain



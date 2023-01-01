#!/bin/bash 

root_dir="."
snps="../data/circadian_variants_introgressed.bed"
iterations=1000
#out_file="../results/temp.txt"

# GENOMATNN
pred="/dors/capra_lab/data/ancient_dna/archaic_hominin/gower21/Nea_to_CEU_af-0.25_w-0.02.tsv"
threshold=0.5
out_file='../data/enrichment_ai_genomatnn.txt'

# MALADAPT
#pred="/dors/capra_lab/data/ancient_dna/archaic_hominin/zhang_MaLAdapt/neanderthal_maladapt_nonafrican.bed"
#threshold=0.9
#out_file='../results/enrichment_ai_maladapt.txt'


python ${root_dir}/enrichment_ai.py $pred $snps $threshold $iterations -o ${out_file}

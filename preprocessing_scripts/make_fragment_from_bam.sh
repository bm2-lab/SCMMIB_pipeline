#!/bin/bash
dir_name=$1

#  single file
sinto fragments -b ${BAM} -p 32 -f site1_donor03_fragments.bed --barcode_regex "[^:]*"


# batch processing
for file_path in $(find "/NFS2_home/NFS2_home_3/wsg/data/BMMC/BAM/$dir_name" -type f \( -name "*.bam" -o -name "*.bam.1" \)); do
    file_name=$(basename "$file_path")
    prefix=$(echo "$file_name" | cut -d'.' -f1)
    nohup /home/wsg/software/miniconda3/envs/sinto/bin/sinto fragments -b $file_path -p 8 -f /NFS2_home/NFS2_home_3/wsg/data/BMMC/Fragments/${dir_name}/${prefix}_fragments.tsv &
done

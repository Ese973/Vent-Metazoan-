# Scripts were ran in Georgia's Advanced Computing Resource Center Sapelo2 Cluster

# Assign Taxonomy 
# Activate QIIME2/2024.10-amplicon
module load QIIME2/2024.10-amplicon

# Script
qiime feature-classifier classify-consensus-blast \
--i-query /directory/deepsea-18s-merged-seqs-05022025.qza \
--i-reference-reads /directory/SILVA138-nr99_fixed_strings_custom_seq_sequences-Mar-10-2025.qza \
--i-reference-taxonomy /directory/SILVA138-nr99_fixed_strings_custom_seq_taxonomy-Mar-10-2025.qza \
--p-maxaccepts 1 \
--p-perc-identity .90 \
--o-classification /directory/18S-rep-sequences-taxonomy.qza \
--o-search-results /directory/analysis/18S-rep-sequences-tax-search-results.qza


# Generating a tree for phylogenetic diversity analyses
# Activate QIIME2/2024.10-amplicon
module load QIIME2/2024.10-amplicon

qiime phylogeny  align-to-tree-mafft-fasttree \
  --i-sequences /directory/denoised-rep-sequences.qza \
  --o-alignment  /directory/aligned-18S-rep-seqs.qza \
  --o-masked-alignment /directory/masked-aligned-rep-seqs.qza \
  --o-tree /directory/unrooted-18S-tree.qza \
  --o-rooted-tree /directory/rooted-18S-tree.qza


# Making nematode phylogeny tree
# Activate QIIME2/2024.10-amplicon
module load QIIME2/2024.10-amplicon

# Script
qiime phylogeny iqtree-ultrafast-bootstrap \
--i-alignment /directory/nema_extra_seqs_aligned_masked.qza \
--p-substitution-model TEST \
--p-alrt 1000 \
--p-n-cores 24 \
--0-tree /directory/nema_extra_bootstrap_tree.qza
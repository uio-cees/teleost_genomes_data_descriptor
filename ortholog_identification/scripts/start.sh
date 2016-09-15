# Michael Matschiner, 2015-03-26

# This pipeline requires the following tools (tested versions in parentheses):
# python3 (Python 3.4.0)
# tblastn (BLAST 2.2.26+)
# mafft (MAFFT v7.157b)
# ruby (Ruby 2.0.0p481)
# java (Java 1.8.0_31)
# perl (Perl v5.18.2)

# Extract the unitigs built from at least 1000 fragments from each of the genome assemblies, these are most likely to represent fragments
# of the mitochondrial genome.
bash reduce_to_mitogenome.sh

# Identification of mitochondrial orthologs.
bash find_orthologs_mitochondrial.sh

# Align dna subject sequences based on amino acid query sequences.
ruby align_dna_subjects.rb

# Filtering of mitochondrial ortholog alignments:
# 1: For each gene, remove poorly aligned sites and sites with too much missing data with BMGE.
ruby filter_sites_with_BMGE.rb

# 2: Remove taxa with too much missing data.
ruby filter_taxa_by_missing_data.rb

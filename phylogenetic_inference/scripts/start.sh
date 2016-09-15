# Michael Matschiner, 2015-03-23.

# This pipeline requires the following tools (tested versions in parentheses):
# raxml (RAxML 8.1.12)
# ruby (Ruby 2.0.0p481)

# Prepare RAxML input files for the concatenated mitochondrial dataset.
ruby prepare_concatenated_RAxML_analyses.rb 

# Run RAxML for the concatenated mitochondrial dataset.
bash run_concatenated_RAxML_analyses.sh

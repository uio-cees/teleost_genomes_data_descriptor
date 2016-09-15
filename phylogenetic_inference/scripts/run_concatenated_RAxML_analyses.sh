# Michael Matschiner, 2015-03-31.

phylogeny_dir="../analysis/mitochondrial/raxml"
for i in ${phylogeny_dir}/*.phy
do
	alignment_file_name=${i}
	partitions_file_name=${i%.phy}.parts
	tree_file_name=${i%.phy}.tre
	python3 resources/run_raxml.py --verbose -b auto -q ${partitions_file_name} ${alignment_file_name} ${tree_file_name}
done
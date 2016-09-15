# Michael Matschiner, 2015-03-26

full_subject_dir="../data/subjects/full"
mitochondrial_subject_dir="../data/subjects/mitochondrial"

for d in ${full_subject_dir}/fish_*.utg.fasta
do
    subject_id=`basename ${d}`
    resources/fastagrep -p "num_frags=\d\d\d\d" $d > ${mitochondrial_subject_dir}/${subject_id}
    makeblastdb -in ${mitochondrial_subject_dir}/${subject_id} -dbtype nucl
done

# Michael Matschiner, 2015-03-26

# Set variables.
query_dir="../data/queries/mitochondrial/protein"
subject_dir="../data/subjects/mitochondrial"
alignment_dir_out="../analysis/alignments/mitochondrial/01"
home_dir=`pwd`

# Make the alignment output directory if it doesn't exist yet.
mkdir -p $alignment_dir_out

# Run find_orthologs.py for each marker, with each of the 76 subjects.
for i in ${query_dir}/*.fasta
do
    marker_id=`basename ${i} | sed 's/.fasta//g'`
    query=`readlink -f ${query_dir}/${marker_id}.fasta`
    mkdir -p ${alignment_dir_out}/${marker_id}
    for i in ${subject_dir}/*.fasta
    do
        subject=`readlink -f ${i}`
        subject_id=`basename $subject | sed 's/.fasta//g' | sed 's/.utg//g'`
        current_out_dir=${alignment_dir_out}/${marker_id}/${subject_id}
        mkdir -p ${current_out_dir}
        cd $current_out_dir
        python3 resources/find_orthologs.py -t -c 2 -e 1e-15 $query $subject
        cd $home_dir	
    done
done

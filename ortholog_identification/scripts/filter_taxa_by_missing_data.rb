## Michael Matschiner 2015-03-20.

# Get the command line arguments.
alignment_directory_in = "../analysis/alignments/mitochondrial/03"
alignment_directory_out = "../analysis/alignments/mitochondrial/04"
minimum_number_of_missing_sequences = 8

# Create the output directory if it does not exist yet.
unless Dir.exists?(alignment_directory_out)
	Dir.mkdir(alignment_directory_out)
end

# Collect names of nucleotide fasta files in the input directory.
mitochondrial_marker_ids = ["ATP6","ATP8","COX1","COX2","COX3","CYTB","ND1","ND2","ND3","ND4","ND4L","ND5"]
dir_entries_in = Dir.entries(alignment_directory_in)
filenames_in = []
dir_entries_in.each do |e|
	if e.match(/.*_nucl.fasta/) or mitochondrial_marker_ids.include?(e.chomp(".fasta"))
		filenames_in << e
	end
end

# Initiate arrays for ids and seqs of all alignments.
ids_per_file = []
seqs_per_file = []

# Do for each fasta file in the input directory.
filenames_in.each do |f|

	# Read the fasta file.
	fasta_file = File.open("#{alignment_directory_in}/#{f}")
	fasta_lines = fasta_file.readlines
	fasta_ids = []
	fasta_bitscores = []
	fasta_hits = []
	fasta_seqs = []
	fasta_lines.each do |l|
		if l[0] == ">"
			fasta_ids << l[1..-1].strip
			fasta_seqs << ""
		else
			fasta_seqs.last << l.strip
		end
	end

	ids_per_file << fasta_ids
	seqs_per_file << fasta_seqs
end

# Get a total list of ids found in all fasta files.
all_ids = []
ids_per_file.each do |i|
	i.each do |ii|
		all_ids << ii unless all_ids.include?(ii)
	end
end
sorted_ids = all_ids.sort

# Count the number of missing sequences per taxon.
missing_per_id = []
sorted_ids.each do |i|
	missing = 0
	ids_per_file.size.times do |x|
		if ids_per_file[x].include?(i)
			missing += 1 if seqs_per_file[x][ids_per_file[x].index(i)].match(/^-+$/)
		else
			missing += 1
		end
	end
	missing_per_id << missing
end

# Get the length of the longes taxon name.
max_id_length = 0
sorted_ids.each do |i|
	max_id_length = i.size if i.size > max_id_length
end

# Report the taxa that are removed from the dataset:
sorted_ids.size.times do |y|
	if missing_per_id[y] > minimum_number_of_missing_sequences
		puts "Excluding #{sorted_ids[y]} with #{missing_per_id[y]} missing sequences."
	end
end

# Write a nexus file for each alignment, only with taxa that have enough data.
ids_per_file.size.times do |x|
	ids_for_this_alignment = []
	seqs_for_this_alignment = []
	sorted_ids.size.times do |y|
		unless missing_per_id[y] > minimum_number_of_missing_sequences
			ids_for_this_alignment << sorted_ids[y]
			if ids_per_file[x].include?(sorted_ids[y])
				seq = seqs_per_file[x][ids_per_file[x].index(sorted_ids[y])]
			else
				seq = ""
				seqs_per_file[x][0].size.times {seq << "-"}
			end
			seqs_for_this_alignment << seq
		end
	end
	nexus_string = "#nexus\n"
	nexus_string << "\n"
	nexus_string << "begin data;\n"
	nexus_string << "dimensions  ntax=#{ids_for_this_alignment.size} nchar=#{seqs_for_this_alignment[0].size};\n"
	nexus_string << "format datatype=DNA gap=- missing=?;\n"
	nexus_string << "matrix\n"
	ids_for_this_alignment.size.times do |y|
		nexus_string << "#{ids_for_this_alignment[y].ljust(max_id_length+2)}#{seqs_for_this_alignment[y]}\n"
	end
	nexus_string << ";\n"
	nexus_string << "end;\n"
	nexus_file_name = "#{filenames_in[x].chomp(".fasta").chomp("_nucl")}.nex"
	nexus_file = File.open("#{alignment_directory_out}/#{nexus_file_name}","w")
	nexus_file.write(nexus_string)
end

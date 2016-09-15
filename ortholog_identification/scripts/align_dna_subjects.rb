## Michael Matschiner 2015-03-26

# Get the command line arguments.
alignment_directory_in = "../analysis/alignments/mitochondrial/01"
alignment_directory_out = "../analysis/alignments/mitochondrial/02"
# The query directory is expected to contain fasta files named [marker_id].fasta
# If such a directory is found, these sequences will be added to the final alignment.
query_directory = "../data/queries/mitochondrial/dna"
# If an exclude list is specified, subject ids in this list will be ignored.
# This list should be formatted like this:
# fish_35:utg7180000000095
# fish_52:utg7180002123755
# ...
# This means that subject sequence with id utg7180000000095 in the subject fish_35.utg.fasta
# (or fish_35.ctg.fasta or fish_35.fasta) will be ignored.
exclude_list_file_name = "../data/subjects/mitochondrial/exclude_list.txt"

# Create the output directory if it does not exist yet.
unless Dir.exists?(alignment_directory_out)
    Dir.mkdir(alignment_directory_out)
end

# Read the exclude list if one has been specified.
exclude_lines = []
if exclude_list_file_name != nil
	exclude_list_file = File.open(exclude_list_file_name)
	exclude_lines = exclude_list_file.readlines
end

# Collect names of marker directories in the input directory.
dir_entries_in = Dir.entries(alignment_directory_in)
marker_dir_names = []
dir_entries_in.each do |e|
	if e[0] != "." 
		unless e.match(".fasta")
			marker_dir_names << e
		end
	end
end

# Go through each marker directory.
marker_dir_names.each do |marker_dir|

	# Collect names of subject directories in the input directory.
	dir_entries_in = Dir.entries("#{alignment_directory_in}/#{marker_dir}")
	subject_dir_names = []
	dir_entries_in.each do |e|
		if e[0] != "."
			unless e.match(".fasta")
				subject_dir_names << e
			end
		end
	end

	# Initiate a string to hold the merged best hits for all subjects for this marker.
	all_best_dna_hits_merged_string = ""

	# Go through each subject directory in the current marker directory.
	subject_dir_names.each do |subject_dir|

		exclude_ids = []
		exclude_lines.each do |l|
			if l.split(":")[0] == subject_dir
				exclude_ids << l.split(":")[1].strip
			end
		end

		# Collect names of alignment files for each query in the subject directory in the current marker directory.
		dir_entries_in = Dir.entries("#{alignment_directory_in}/#{marker_dir}/#{subject_dir}")
		query_file_names = []
		dir_entries_in.each do |e|
			if e.match(/.fasta/) and e[0..2] != "all" and e[0..3] != "best"
				query_file_names << e
			end
		end

		# Sort fasta file names into dna and protein fasta files.
		dna_query_file_names = []
		protein_query_file_names = []
		query_file_names.each do |f|
			if f.match(/_nucl.fasta/)
				dna_query_file_names << f
			else
				protein_query_file_names << f
			end
		end

		# Read all protein query files and collect the first sequence from each of these.
		ids = []
		seqs = []
		protein_query_file_names.each do |q|
			protein_query_file = File.open("#{alignment_directory_in}/#{marker_dir}/#{subject_dir}/#{q}")
			protein_query_lines = protein_query_file.readlines
			protein_query_file.close
			ids << protein_query_lines[0][1..-1].strip
			seqs << protein_query_lines[1].strip
		end

		# Prepare a new fasta string.
		all_protein_queries_string = ""
		ids.size.times do |x|
			all_protein_queries_string << ">#{ids[x]}\n"
			# Replace existing gaps with 'X' as these should not be removed by mafft anyway.
			all_protein_queries_string << "#{seqs[x].gsub("-","X")}\n"
		end

		# Write a new fasta file containing only the first sequence of each protein query file
		# (the other sequences contain hits that were obtained with this query).
		all_protein_queries_file_name = "#{alignment_directory_in}/#{marker_dir}/#{subject_dir}/all_protein_queries.fasta"
		all_protein_queries_aligned_file_name =  "#{all_protein_queries_file_name.chomp(".fasta")}_aln.fasta"
		all_protein_queries_file = File.open(all_protein_queries_file_name,"w")
		all_protein_queries_file.write(all_protein_queries_string)
		all_protein_queries_file.close

		# Run mafft to align the fasta file containing all protein queries for this combination of subject and marker.
		system("mafft #{all_protein_queries_file_name} > #{all_protein_queries_aligned_file_name} 2> /dev/null")

		# Read the fasta file aligned by mafft.
		all_protein_queries_aligned_file = File.open(all_protein_queries_aligned_file_name)
		all_protein_queries_aligned_lines = all_protein_queries_aligned_file.readlines
		all_protein_queries_aligned_file.close
		ids = []
		protein_seqs = []
		gaps_per_protein_seq = []
		all_protein_queries_aligned_lines.each do |l|
			if l[0] == ">"
				ids << l[1..-1].strip
				protein_seqs << ""
			elsif l.strip != ""
				protein_seqs.last << l.strip
			end
		end

		# Produce arrays for the gaps contained in each sequence.
		protein_seqs.each do |s|
			gaps = []
			s.size.times do |x|
				gaps << x if s[x] == "-"
			end
			gaps_per_protein_seq << gaps
		end

		# Read all hits for the protein queries for this combination of subject and marker.
		ids_per_protein_query_file = []
		protein_seqs_per_protein_query_file = []
		protein_query_file_names.size.times do |x|
			ids = []
			protein_seqs = []
			protein_query_file = File.open("#{alignment_directory_in}/#{marker_dir}/#{subject_dir}/#{protein_query_file_names[x]}")
			protein_query_lines = protein_query_file.readlines
			protein_query_file.close
			protein_query_lines[2..-1].each do |l|
				if l[0] == ">"
					unless l[1..-1].include?("sseqid=None")
						ids << l[1..-1].strip
						protein_seqs << ""
					end
				elsif l.strip != "" and l.strip != "-"
					protein_seqs.last << l.strip
				end
			end
			ids.size.times do |y|
				ids[y].match(/sseqid=(.+?)[,\]]/)
				sseqid=$1
				if exclude_ids.include?(sseqid)
					ids[y] = nil
					protein_seqs[y] = nil
				end
			end
			ids_per_protein_query_file << ids.compact
			protein_seqs_per_protein_query_file << protein_seqs.compact
		end

		# Read all hits for the dna queries for this combination of subject and marker.
		ids_per_dna_query_file = []
		dna_seqs_per_dna_query_file = []
		dna_query_file_names.size.times do |x|
			ids = []
			dna_seqs = []
			dna_query_file = File.open("#{alignment_directory_in}/#{marker_dir}/#{subject_dir}/#{dna_query_file_names[x]}")
			dna_query_lines = dna_query_file.readlines
			dna_query_file.close
			dna_query_lines.each do |l|
				if l[0] == ">"
					ids << l[1..-1].strip
					dna_seqs << ""
				elsif l.strip != ""
					dna_seqs.last << l.strip
				end
			end
			ids.size.times do |y|
				ids[y].match(/sseqid=(.+?)[,\]]/)
				sseqid=$1
				if exclude_ids.include?(sseqid)
					ids[y] = nil
					dna_seqs[y] = nil
				end
			end
			ids_per_dna_query_file << ids.compact
			dna_seqs_per_dna_query_file << dna_seqs.compact
		end

		# Make sure that the total number of protein sequences matches that of dna sequences.
		total_number_of_protein_sequences = 0
		total_number_of_dna_sequences = 0
		unless protein_seqs_per_protein_query_file.size == dna_seqs_per_dna_query_file.size
			raise "ERROR: The numbers of protein and dna sequence files don't match for subject #{subject_dir}!"
		end
		protein_seqs_per_protein_query_file.size.times do |x|
			if dna_seqs_per_dna_query_file[x].size != protein_seqs_per_protein_query_file[x].size
				raise "ERROR: The number of dna and protein sequences do not match for #{ids_per_dna_query_file[x]} (x = #{x}) and subject #{subject_dir}!"
			end
		end

		# Make sure that each dna sequence is three times the length of the protein sequence.
		protein_seqs_per_protein_query_file.size.times do |x|
			protein_seqs_per_protein_query_file[x].size.times do |y|
				if protein_seqs_per_protein_query_file[x][y].size * 3 != dna_seqs_per_dna_query_file[x][y].size
					# If the length don't match, check whether the lengths would check after removing all gaps.
					if protein_seqs_per_protein_query_file[x][y].gsub("-","").size * 3 == dna_seqs_per_dna_query_file[x][y].gsub("-","").size
						# Fix the dna sequence by adding gaps.
						old_dna_seq = dna_seqs_per_dna_query_file[x][y].gsub("-","")
						new_dna_seq = ""
						protein_seqs_per_protein_query_file[x][y].size.times do |pos|
							if protein_seqs_per_protein_query_file[x][y][pos] == "-"
								new_dna_seq << "---"
							else
								new_dna_seq << old_dna_seq.slice!(0..2)
							end
						end
						dna_seqs_per_dna_query_file[x][y] = new_dna_seq
					else
						puts "ERROR: The length of dna sequence #{ids_per_dna_query_file[x][y]} does not match that of the respective protein sequence for file #{protein_query_file_names[x]} and subject #{subject_dir}!"
						puts "dna:      #{dna_seqs_per_dna_query_file[x][y]}"
						puts "protein:  #{protein_seqs_per_protein_query_file[x][y]}"
						raise
					end
				end
			end
		end

		# Use the array of gaps produced above to align the protein and dna hit sequences.
		aligned_protein_seqs_per_protein_query_file = []
		aligned_dna_seqs_per_protein_query_file = []
		protein_query_file_names.size.times do |x|
			aligned_protein_seqs = []
			aligned_dna_seqs = []
			protein_seqs_per_protein_query_file[x].size.times do |y|
				aligned_protein_seq = ""
				aligned_dna_seq = ""
				pos_in_seq = 0
				pos_in_alignment = 0
				while pos_in_seq <= protein_seqs_per_protein_query_file[x][y].size
					if gaps_per_protein_seq[x].include?(pos_in_alignment)
						aligned_protein_seq << "-"
						aligned_dna_seq << "---"
						pos_in_alignment += 1
					else
						if pos_in_seq < protein_seqs_per_protein_query_file[x][y].size
							aligned_protein_seq << protein_seqs_per_protein_query_file[x][y][pos_in_seq]
							aligned_dna_seq << dna_seqs_per_dna_query_file[x][y][pos_in_seq*3..pos_in_seq*3+2]
						end
						pos_in_seq += 1
						pos_in_alignment += 1
					end
				end
				aligned_protein_seqs << aligned_protein_seq
				aligned_dna_seqs << aligned_dna_seq
			end
			aligned_protein_seqs_per_protein_query_file << aligned_protein_seqs
			aligned_dna_seqs_per_protein_query_file << aligned_dna_seqs
		end

		# Prepare a string for the alignment of protein hits for this combination of subject and marker.
		all_protein_hits_string = ""
		ids_per_protein_query_file.size.times do |x|
			ids_per_protein_query_file[x].size.times do |y|
				all_protein_hits_string << ">#{ids_per_protein_query_file[x][y]}\n"
				all_protein_hits_string << "#{aligned_protein_seqs_per_protein_query_file[x][y]}\n"
			end
		end

		# Prepare a string for the alignment of dna hits for this combination of subject and marker.
		all_dna_hits_string = ""
		ids_per_dna_query_file.size.times do |x|
			ids_per_dna_query_file[x].size.times do |y|
				all_dna_hits_string << ">#{ids_per_dna_query_file[x][y]}\n"
				all_dna_hits_string << "#{aligned_dna_seqs_per_protein_query_file[x][y]}\n"
			end
		end		

		# Write the string for the alignment of protein hits for this combination of subject and marker to a fasta file.
		all_protein_hits_file_name = "#{alignment_directory_in}/#{marker_dir}/#{subject_dir}/all_protein_hits.fasta"
		all_protein_hits_file = File.open(all_protein_hits_file_name,"w")
		all_protein_hits_file.write(all_protein_hits_string)
		all_protein_hits_file.close

		# Write the string for the alignment of dna hits for this combination of subject and marker to a fasta file.
		all_dna_hits_file_name = "#{alignment_directory_in}/#{marker_dir}/#{subject_dir}/all_dna_hits.fasta"
		all_dna_hits_file = File.open(all_dna_hits_file_name,"w")
		all_dna_hits_file.write(all_dna_hits_string)
		all_dna_hits_file.close

		# For each unitig id, calculate their total (=sum of all) bitscore.
		hit_unitig_ids = []
		hit_unitig_total_bitscores = []
		ids_per_dna_query_file.size.times do |x|
			ids_per_dna_query_file[x].size.times do |y|
				# Get the unitig id of this hit.
				ids_per_dna_query_file[x][y].match(/&sseqid=(utg\d+),/)
				hit_unitig_id = $1
				if hit_unitig_id == nil
					ids_per_dna_query_file[x][y].match(/&sseqid=(\w{6}),/)
					hit_unitig_id = $1
				end
				# Store this unitig id if it isn't already.
				unless hit_unitig_ids.include?(hit_unitig_id)
					hit_unitig_ids << hit_unitig_id
					hit_unitig_total_bitscores << 0
				end
				# Get the bitscore of this hit.
				ids_per_dna_query_file[x][y].match(/bitscore=(\d+)\./)
				hit_unitig_bitscore = $1.to_i
				# Add this hits' bitscore to the respective unitig id's bitscore.
				hit_unitig_total_bitscores[hit_unitig_ids.index(hit_unitig_id)] += hit_unitig_bitscore
			end
		end

		# For each unitig id, find the hit with the best bitscore.
		hit_unitig_best_ids = []
		hit_unitig_best_seqs = []
		hit_unitig_best_bitscores = []
		hit_unitig_ids.each do |selected_hit_unitig_id|
			hit_unitig_best_id = ""
			hit_unitig_best_seq = ""
			hit_unitig_best_bitscore = 0
			ids_per_dna_query_file.size.times do |x|
				ids_per_dna_query_file[x].size.times do |y|
					# Get the unitig id of this hit.
					ids_per_dna_query_file[x][y].match(/&sseqid=(utg\d+),/)
					hit_unitig_id = $1
					if hit_unitig_id == nil
						ids_per_dna_query_file[x][y].match(/&sseqid=(\w{6}),/)
						hit_unitig_id = $1
					end
					# If the id of this unitig matches the currently selected unitig id.
					if hit_unitig_id == selected_hit_unitig_id
						# Get the bitscore of this hit.
						ids_per_dna_query_file[x][y].match(/bitscore=(\d+)\./)
						hit_unitig_bitscore = $1.to_i
						# If this unitig's bitscore is better than all previous ones,
						# store its id and sequence.
						if hit_unitig_bitscore > hit_unitig_best_bitscore
							hit_unitig_best_id = ids_per_dna_query_file[x][y]
							hit_unitig_best_seq = aligned_dna_seqs_per_protein_query_file[x][y]
							hit_unitig_best_bitscore = hit_unitig_bitscore
						end
					end
				end
			end
			hit_unitig_best_ids << hit_unitig_best_id
			hit_unitig_best_seqs << hit_unitig_best_seq
			hit_unitig_best_bitscores << hit_unitig_best_bitscore
		end

		# Sort the arrays ids and sequences of hits with the best bitscores by their bitscores.
		sorted = false
		until sorted
			sorted = true
			0.upto(hit_unitig_best_ids.size-2) do |x|
				(x+1).upto(hit_unitig_best_ids.size-1) do |y|
					if hit_unitig_best_bitscores[y] > hit_unitig_best_bitscores[x]
						hit_unitig_best_ids[x],hit_unitig_best_ids[y] = hit_unitig_best_ids[y],hit_unitig_best_ids[x]
						hit_unitig_best_seqs[x],hit_unitig_best_seqs[y] = hit_unitig_best_seqs[y],hit_unitig_best_seqs[x]
						hit_unitig_best_bitscores[x],hit_unitig_best_bitscores[y] = hit_unitig_best_bitscores[y],hit_unitig_best_bitscores[x]
						sorted = false
					end
				end
			end
		end

		# Merge hits of different unitigs if their sequences overlap on at least 100 bp and are identical within the overlapping region.
		merged_completely = false
		hit_unitig_best_ids_merged = Marshal.load(Marshal.dump(hit_unitig_best_ids))
		hit_unitig_best_seqs_merged = Marshal.load(Marshal.dump(hit_unitig_best_seqs))
		until merged_completely
			merged_completely = true
			0.upto(hit_unitig_best_ids_merged.size-2) do |x|
				if merged_completely
					(x+1).upto(hit_unitig_best_ids_merged.size-1) do |y|
						if merged_completely

							# Test whether the two sequence overlap and are identical within the overlap.
							overlap_length = 0
							overlap_identity = true
							hit_unitig_best_seqs_merged[x].size.times do |pos|
								if hit_unitig_best_seqs_merged[x][pos] != "-" and hit_unitig_best_seqs_merged[y][pos] != "-"
									if hit_unitig_best_seqs_merged[y][pos] == hit_unitig_best_seqs_merged[x][pos]
										overlap_length += 1
									else
										overlap_identity = false
									end
								end
							end
							# If the overlap is identical and sufficiently long, merge the two sequences.
							if overlap_length >= 100 and overlap_identity == true
								id_string1 = hit_unitig_best_ids_merged[x]
								id_string2 = hit_unitig_best_ids_merged[y]
								id1 = id_string1.sub(/\[.+\]/,"")
								id2 = id_string2.sub(/\[.+\]/,"")
								if id1 == id2
									id_merged = id1
								else
									raise "ERROR: Found different ids in sequences to be merged!"
								end
								id_string1.match(/sseqid=(.+?)[,\]]/)
								sseqid1 = $1
								id_string2.match(/sseqid=(.+?)[,\]]/)
								sseqid2 = $1
								id_string1.match(/bitscore=(.+?)[,\]]/)
								bitscore1 = $1.to_f
								id_string2.match(/bitscore=(.+?)[,\]]/)
								bitscore2 = $1.to_f
								sequence_length1 = hit_unitig_best_seqs_merged[x].gsub("-","").size
								sequence_length2 = hit_unitig_best_seqs_merged[y].gsub("-","").size
								if bitscore1 > bitscore2
									bitscore_merged = (bitscore1 + bitscore2 * ((sequence_length2-overlap_length)/sequence_length2.to_f)).round
								else
									bitscore_merged = (bitscore2 + bitscore1 * ((sequence_length1-overlap_length)/sequence_length1.to_f)).round
								end
								sseqid_merged = "#{sseqid1}_#{sseqid2}"
								id_string_merged = "#{id_merged}[&sseqid=#{sseqid_merged},bitscore=#{bitscore_merged}.0]"
								seq_merged = ""
								hit_unitig_best_seqs_merged[x].size.times do |pos|
									if hit_unitig_best_seqs_merged[x][pos] != "-"
										seq_merged << hit_unitig_best_seqs_merged[x][pos]
									else
										seq_merged << hit_unitig_best_seqs_merged[y][pos]
									end
								end
								tmp_ids = []
								tmp_seqs = []
								hit_unitig_best_ids_merged.size.times do |z|
									if z == x
										tmp_ids << id_string_merged
										tmp_seqs << seq_merged
									elsif z != y
										tmp_ids << hit_unitig_best_ids_merged[z]
										tmp_seqs << hit_unitig_best_seqs_merged[z]
									end
								end
								hit_unitig_best_ids_merged = Marshal.load(Marshal.dump(tmp_ids))
								hit_unitig_best_seqs_merged = Marshal.load(Marshal.dump(tmp_seqs))
								merged_completely = false
							end
						end
					end
				end
			end
		end

		# Prepare a string for the alignment of the best dna hits per unitig for this combination of subject and marker.
		best_dna_hits_string = ""
		hit_unitig_best_ids.size.times do |x|
			best_dna_hits_string << ">#{hit_unitig_best_ids[x]}\n"
			best_dna_hits_string << "#{hit_unitig_best_seqs[x]}\n"
		end		

		# Write the string for the alignment of the best dna hits per unitig for this combination of subject and marker to a fasta file.
		best_dna_hits_file_name = "#{alignment_directory_in}/#{marker_dir}/#{subject_dir}/best_dna_hits.fasta"
		best_dna_hits_file = File.open(best_dna_hits_file_name,"w")
		best_dna_hits_file.write(best_dna_hits_string)
		best_dna_hits_file.close

		# Prepare a string for the alignment of the merged best dna hits per unitig for this combination of subject and marker.
		best_dna_hits_merged_string = ""
		hit_unitig_best_ids_merged.size.times do |x|
			best_dna_hits_merged_string << ">#{hit_unitig_best_ids_merged[x]}\n"
			all_best_dna_hits_merged_string << ">#{hit_unitig_best_ids_merged[x]}\n"
			best_dna_hits_merged_string << "#{hit_unitig_best_seqs_merged[x].downcase}\n"
			all_best_dna_hits_merged_string << "#{hit_unitig_best_seqs_merged[x].downcase.gsub("-","")}\n"
		end		

		# Write the string for the alignment of the best merged dna hits per unitig for this combination of subject and marker to a fasta file.
		best_dna_hits_merged_file_name = "#{alignment_directory_in}/#{marker_dir}/#{subject_dir}/best_dna_hits_merged.fasta"
		best_dna_hits_merged_file = File.open(best_dna_hits_merged_file_name,"w")
		best_dna_hits_merged_file.write(best_dna_hits_merged_string)
		best_dna_hits_merged_file.close

	end

	# Write the string for the alignment of all best merged dna hits per unitig for this marker to a fasta file.
	all_best_dna_hits_merged_file_name = "#{alignment_directory_in}/#{marker_dir}/all_best_dna_hits_merged.fasta"
	all_best_dna_hits_merged_file = File.open(all_best_dna_hits_merged_file_name,"w")
	all_best_dna_hits_merged_file.write(all_best_dna_hits_merged_string)
	all_best_dna_hits_merged_file.close

	# Combine this fasta file with the fasta file containing all query sequences for this marker.
	original_query_file_name = "#{query_directory}/#{marker_dir}.fasta"
	all_best_dna_hits_merged_and_queries_file_name = "#{alignment_directory_in}/#{marker_dir}/all_best_dna_hits_merged_and_queries.fasta"
	system("cat #{all_best_dna_hits_merged_file_name} #{original_query_file_name} > #{all_best_dna_hits_merged_and_queries_file_name}")

	# Run translatorx to align the fasta file containing all protein queries for this combination of subject and marker.
	system("perl resources/translatorx.pl -i #{all_best_dna_hits_merged_and_queries_file_name} -c 2 -p F 2> /dev/null")
	# system("mafft #{all_best_dna_hits_merged_and_queries_file_name} > #{all_best_dna_hits_merged_and_queries_aln_file_name} 2> /dev/null")
	
	# Read the translatorx output dna alignment.
	# all_best_dna_hits_merged_and_queries_aln_file_name = "#{all_best_dna_hits_merged_and_queries_file_name.chomp(".fasta")}_aln.fasta"
	translatorx_res_file = File.open("translatorx_res.nt_ali.fasta")
	translatorx_res_lines = translatorx_res_file.readlines
	translatorx_res_file.close

	# Delete all output from translatorx.
	system("rm translatorx_res.*")

	# Analyse the translatorx output dna alignment.
	translatorx_res_ids = []
	translatorx_res_seqs = []
	translatorx_res_lines.each do |l|
		if l[0] == ">"
			translatorx_res_id = l[1..-1].strip
			if translatorx_res_id.include?("&")
				unless translatorx_res_id.include?("[&")
					translatorx_res_id.sub!("&","[&")
					translatorx_res_id << "]"
				end
			end
			translatorx_res_ids << translatorx_res_id
			translatorx_res_seqs << ""
		elsif l.strip != ""
			translatorx_res_seqs.last << l.strip
		end
	end

	# Prepare the string for the alignment of all best merged dna hits per unitig plus queries for this marker.
	all_best_dna_hits_merged_and_queries_string = ""
	translatorx_res_ids.size.times do |x|
		all_best_dna_hits_merged_and_queries_string << ">#{translatorx_res_ids[x]}\n"
		all_best_dna_hits_merged_and_queries_string << "#{translatorx_res_seqs[x]}\n"
	end

	# Write the string for the alignment of all best merged dna hits per unitig plus queries for this marker to a fasta file.
	all_best_dna_hits_merged_and_queries_file_name = "#{alignment_directory_in}/#{marker_dir}/all_best_dna_hits_merged_and_queries_aln.fasta"
	all_best_dna_hits_merged_and_queries_file = File.open(all_best_dna_hits_merged_and_queries_file_name,"w")
	all_best_dna_hits_merged_and_queries_file.write(all_best_dna_hits_merged_and_queries_string)
	all_best_dna_hits_merged_and_queries_file.close

	# Get a list of all taxon ids.
	all_taxon_ids = []
	translatorx_res_ids.each do |i|
		taxon_id = i.sub(/\[.*\]/,"")
		all_taxon_ids << taxon_id unless all_taxon_ids.include?(taxon_id)
	end
	sorted_taxon_ids = all_taxon_ids.sort

	# For each taxon id, get all sequences of it and build a consensus.
	consensus_seqs = []
	sorted_taxon_ids.each do |i|
		taxon_seqs = []
		translatorx_res_ids.size.times do |x|
			if translatorx_res_ids[x].sub(/\[.*\]/,"") == i
				taxon_seqs << translatorx_res_seqs[x]
			end
		end
		taxon_consensus = ""
		taxon_seqs[0].size.times do |pos|
			bases_at_this_pos = []
			taxon_seqs.each do |s|
				bases_at_this_pos << s[pos]
			end
			bases_at_this_pos.uniq!
			bases_at_this_pos.delete("-")
			if bases_at_this_pos.size == 0
				taxon_consensus << "-"
			elsif bases_at_this_pos.size == 1
				taxon_consensus << bases_at_this_pos[0]
			else
				taxon_consensus << "n"
			end
		end
		consensus_seqs << taxon_consensus
	end

	# Prepare a string for the output alignment with all consensus sequences.
	consensus_string = ""
	sorted_taxon_ids.size.times do |x|
		consensus_string << ">#{sorted_taxon_ids[x]}\n"
		consensus_string << "#{consensus_seqs[x]}\n"
	end

	# Write the consensus sequence string to a file in the alignment output directory.
	consensus_file_name = "#{alignment_directory_out}/#{marker_dir}.fasta"
	consensus_file = File.open(consensus_file_name,"w")
	consensus_file.write(consensus_string)
	consensus_file.close

end


function get_seq_range(record::FASTQ.Record)::UnitRange{Int64}
	(record.description_len + 1):(record.description_len + Int64(record.has_description_seq_len))
end


function get_qual_range(record::FASTQ.Record)::UnitRange{Int64}
	(record.description_len + Int64(record.has_description_seq_len) + 1):length(record.data)
end


get_seq(record::FASTQ.Record)::Vector{UInt8} = record.data[get_seq_range(record)]


get_rev_comp(record::FASTQ.Record)::Vector{UInt8} = reverse(replace(record.data[get_seq_range(record)], 0x54 => 0x41, 0x41 => 0x54, 0x43 => 0x47, 0x47 => 0x43))


get_qual(record::FASTQ.Record)::Vector{UInt8} = record.data[get_qual_range(record)]


function sliding_window_kmer_indices(seq::Vector{UInt8}, k::Int)
    kmer_dict = Dict{Vector{UInt8}, Vector{Int}}()
    for i in 1:(length(seq) - k + 1)
        kmer = seq[i:i+k-1]
        push!(get!(kmer_dict, kmer, []), i)
    end
    return Dict(key => val[1] for (key, val) in kmer_dict if length(val) == 1)
end


function get_kmer_offsets(forward::FASTQ.Record, reverse::FASTQ.Record, window_size::Int64)::Int64
	r1_kmers = sliding_window_kmer_indices(get_seq(forward), window_size)
	r2_kmers = sliding_window_kmer_indices(get_rev_comp(reverse), window_size)
	common_kmers = intersect(keys(r1_kmers), keys(r2_kmers))
	count_map = countmap(r1_kmers[kmer] - r2_kmers[kmer] for kmer in common_kmers)
	total_count = sum(values(count_map))
	relative_map = Dict(k => v / total_count for (k, v) in count_map)
	if maximum(values(relative_map)) < 0.5
		error("Unable to find an unambiguous alignment.")
	end
	findmax(relative_map)[2]
end


function align_sequences(forw::FASTQ.Record, rev::FASTQ.Record)::Matrix{UInt8}
	for_seq = get_seq(forw)
	rev_seq = get_rev_comp(rev)
	for_qual = get_qual(forw)
	rev_qual = reverse(get_qual(rev))
	for_len = length(for_seq)
	rev_len = length(rev_seq)
	offset = get_kmer_offsets(forw, rev, 5)
	if offset <= 0
		max_len = maximum((rev_len, for_len - offset))
		for_ix = 1 - offset:for_len - offset 
		rev_ix = 1:rev_len
	elseif offset > 0
		max_len = maximum((for_len, rev_len + offset))
		for_ix = 1:for_len
		rev_ix = 1 + offset:rev_len + offset
	end
	alignment_matrix = fill(0x21, max_len, 4)
	alignment_matrix[for_ix, 1] .= for_seq
	alignment_matrix[rev_ix, 2] .= rev_seq
	alignment_matrix[for_ix, 3] .= for_qual
	alignment_matrix[rev_ix, 4] .= rev_qual
	return alignment_matrix
end
    

function calculate_ds_qual(read_matr::Matrix{UInt8})
	read_matr = Float32.(read_matr .- 0x21)
	read_matr[:, 3] .= 10 .^ -(read_matr[:, 3] ./ 10)
	read_matr[:, 4] .= 10 .^ -(read_matr[:, 4] ./ 10)
	out_matr = zeros(Float32, size(read_matr, 1), 2)
	for i in 1:size(read_matr, 1)
		out_matr[i, :] .= calculate_posteriori(read_matr[i, :]...)
	end
	out_matr[:, 2] .= round.(-10 .* log10.(out_matr[:, 2]))
	out_matr = UInt8.(clamp.(out_matr, 0, 255) .+ 0x21)
	return out_matr
end


function calculate_posteriori(for_base, rev_base, for_prob, rev_prob)
	if for_prob == 1.0
		posteriori = rev_prob
		ret_base = rev_base
	elseif rev_prob == 1.0
		posteriori = for_prob
		ret_base = for_base
	elseif for_base == rev_base
		posteriori = ((for_prob * rev_prob / 3.0)
			/ (1.0 - for_prob - rev_prob + 4.0
				* for_prob * rev_prob / 3.0))  # (8) in paper
		ret_base = for_base
	else
		if for_prob < rev_prob
			posteriori = ((for_prob * (1.0 - rev_prob / 3.0)
				/ (for_prob + rev_prob - 4.0
					* for_prob * rev_prob / 3.0)))  # (9) in paper
			ret_base = for_base
		else
			posteriori = ((rev_prob * (1.0 - for_prob / 3.0)
				/ (for_prob + rev_prob - 4.0
					* for_prob * rev_prob / 3.0)))  # (9) in paper
			ret_base = rev_base
		end
	end
	return ret_base, posteriori
end


function merge_fastq_reads(for_read::FASTX.FASTQ.Record, rev_read::FASTX.FASTQ.Record)::FASTX.FASTQ.Record
	id_string = identifier(for_read)
	aligned = align_sequences(for_read, rev_read)
	merged = calculate_ds_qual(aligned)
	return FASTX.FASTQ.Record(id_string, String(merged[:,1]), String(merged[:,2]))
end


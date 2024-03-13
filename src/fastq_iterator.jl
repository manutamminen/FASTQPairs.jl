
mutable struct PairedFastqIterator
	reader1::FASTQ.Reader
	reader2::FASTQ.Reader
	ix::Int64
end


function PairedFastqIterator(forward_reads::String, reverse_reads::String)
    stream1 = endswith(forward_reads, ".gz") ? GzipDecompressorStream(open(forward_reads)) : open(forward_reads)
    stream2 = endswith(reverse_reads, ".gz") ? GzipDecompressorStream(open(reverse_reads)) : open(reverse_reads)
    return PairedFastqIterator(FASTQ.Reader(stream1), FASTQ.Reader(stream2), 1)
end


function Base.iterate(iter::PairedFastqIterator, state=nothing)
	read1 = iterate(iter.reader1, state)
	read2 = iterate(iter.reader2, state)
	if read1 === nothing && read2 === nothing
		return nothing
	elseif read1 === nothing && read2 !== nothing
		close(iter.reader2)
		error("Reverse read file longer than forward read file")
	elseif read1 !== nothing && read2 === nothing
		close(iter.reader1)
		error("Forward read file longer than reverse read file")
	end
	id_vect = read1[1].data[1:read1[1].description_len] - read2[1].data[1:read2[1].description_len]
	if sum(id_vect .!= 0) != 1
		close(iter.reader1)
		close(iter.reader2)
		error("Mismatching Fastq ids at $(iter.ix)!")
	end
	iter.ix += 1
	return (read1[1], read2[1]), (read1[2], read2[2])
end


Base.IteratorSize(::PairedFastqIterator) = Base.SizeUnknown()


Base.eltype(::PairedFastqIterator) = Tuple{FASTX.FASTQ.Record, FASTX.FASTQ.Record}


# struct BatchPairedFastqIterator
#     iterator::PairedFastqIterator
#     batch_size::Int
# end


# function Base.iterate(iter::BatchPairedFastqIterator, state=(nothing, 1))
#     next_state, batch_index = state
#     batch = []
#     count = 0
#     while count < iter.batch_size
#         item = iterate(iter.iterator, next_state)
#         if item === nothing
#             break
#         end
#         value, next_state = item
#         push!(batch, value)
#         count += 1
#     end
#     if isempty(batch)
#         return nothing
#     end
#     return (batch, (next_state, batch_index + 1))
# end


# function BatchPairedFastqIterator(forward_reads::String, reverse_reads::String, batch_size::Int)
# 	paired_fq_iter = PairedFastqIterator(forward_reads, reverse_reads)
# 	return BatchPairedFastqIterator(paired_fq_iter, batch_size)
# end


struct BatchPairedFastqIterator
    iterator::PairedFastqIterator
    batch_size::Int
    buffer::Vector{Tuple{FASTX.FASTQ.Record, FASTX.FASTQ.Record}}  # Pre-allocated buffer
end


function BatchPairedFastqIterator(forward_reads::String, reverse_reads::String, batch_size::Int)
    paired_fq_iter = PairedFastqIterator(forward_reads, reverse_reads)
    buffer = Vector{Tuple{FASTX.FASTQ.Record, FASTX.FASTQ.Record}}(undef, batch_size)
    return BatchPairedFastqIterator(paired_fq_iter, batch_size, buffer)
end


function Base.iterate(iter::BatchPairedFastqIterator, state=(nothing, 1))
    next_state, batch_index = state
    count = 0
    while count < iter.batch_size
        item = iterate(iter.iterator, next_state)
        if item === nothing
            break
        end
        value, next_state = item
        iter.buffer[count + 1] = value
        count += 1
    end
    if count == 0
        return nothing
    end
    return (view(iter.buffer, 1:count), (next_state, batch_index + 1))
end

Base.IteratorSize(::BatchPairedFastqIterator) = Base.SizeUnknown()


Base.eltype(::BatchPairedFastqIterator) = Vector{Tuple{FASTX.FASTQ.Record, FASTX.FASTQ.Record}}


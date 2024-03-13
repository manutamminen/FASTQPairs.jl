module FASTQPairs

using FASTX, BioSequences, CodecZlib, StatsBase

export PairedFastqIterator, BatchPairedFastqIterator, merge_fastq_reads


include("fastq_iterator.jl")
include("read_merge.jl")


end

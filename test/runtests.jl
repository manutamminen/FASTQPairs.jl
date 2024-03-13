using Revise
using FASTQPairs
using Test


@testset "FASTQPairs.jl" begin
    corrects = collect(PairedFastqIterator("data/correct_R1.fastq.gz", "data/correct_R2.fastq.gz"))
    @test length(corrects) == 8
    uncompressed = collect(PairedFastqIterator("data/uncompressed_R1.fastq", "data/uncompressed_R2.fastq"))
    @test length(uncompressed) == 3
    incorrect_for = PairedFastqIterator("data/incorrect_R1.fastq.gz", "data/correct_R2.fastq.gz")
    incorrect_rev = PairedFastqIterator("data/correct_R1.fastq.gz", "data/incorrect_R2.fastq.gz")
    long_for = PairedFastqIterator("data/long_correct_R1.fastq.gz", "data/correct_R2.fastq.gz")
    long_rev = PairedFastqIterator("data/correct_R1.fastq.gz", "data/long_correct_R2.fastq.gz")
    twice_for = PairedFastqIterator("data/correct_R1.fastq.gz", "data/correct_R1.fastq.gz")
    @test_throws ErrorException collect(incorrect_for)
    @test_throws ErrorException collect(incorrect_rev) 
    @test_throws ErrorException collect(long_for) 
    @test_throws ErrorException collect(long_rev) 
    @test_throws ErrorException collect(twice_for) 
    @test FASTQPairs.get_kmer_offsets(corrects[1]..., 5) == 0
    @test FASTQPairs.get_kmer_offsets(corrects[2]..., 5) == 1
    @test FASTQPairs.get_kmer_offsets(corrects[3]..., 5) == -1
    @test FASTQPairs.get_kmer_offsets(corrects[5]..., 5) == -84
    @test FASTQPairs.get_kmer_offsets(corrects[6]..., 5) == -85
    @test FASTQPairs.get_kmer_offsets(corrects[7]..., 5) == -83
    @test_throws ErrorException FASTQPairs.get_kmer_offsets(corrects[4]..., 5)
    @test_throws ErrorException FASTQPairs.get_kmer_offsets(corrects[8]..., 5)
    @test size(FASTQPairs.align_sequences(corrects[1]...)) == (151, 4)
    @test size(FASTQPairs.align_sequences(corrects[2]...)) == (152, 4)
    @test size(FASTQPairs.align_sequences(corrects[3]...)) == (151, 4)
    @test size(FASTQPairs.align_sequences(corrects[5]...)) == (151, 4)
    @test size(FASTQPairs.align_sequences(corrects[6]...)) == (152, 4)
    @test size(FASTQPairs.align_sequences(corrects[7]...)) == (151, 4)
    @test_throws ErrorException FASTQPairs.align_sequences(corrects[4]...)
    @test_throws ErrorException FASTQPairs.align_sequences(corrects[8]...)
end


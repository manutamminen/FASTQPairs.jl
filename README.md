# FASTQPairs.jl

[![Build Status](https://github.com/manutamminen/FASTQPairs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/manutamminen/FASTQPairs.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Introduction
FASTQPairs.jl is a Julia package developed for efficient processing of paired-end FASTQ files. Designed with bioinformatics workflows in mind, this tool offers a streamlined approach for handling, analyzing, and manipulating high-throughput sequencing data in the FASTQ format. Whether you're working on genome assembly, variant calling, or transcriptome analysis, FASTQPairs.jl provides the necessary functionality to work effectively with paired-end sequencing data.

## Features
- Easy-to-use functions for reading and processing paired FASTQ files
- Support for both synchronous and asynchronous processing to optimize performance
- Integration capabilities with other Julia bioinformatics tools for comprehensive analysis pipelines

## Installation
You can install FASTQPairs.jl from the julia REPL. Press ] to enter pkg mode, and enter the following:

```julia
add https://github.com/manutamminen/FASTQPairs.jl.git
```

## Usage
After installation, you can start using FASTQPairs.jl in your Julia scripts or interactive sessions. Here's a quick example to get you started:

```julia
using FASTQPairs

sequence_iterator = PairedFastqIterator("data/correct_R1.fastq.gz", "data/correct_R2.fastq.gz")

for sequence_pair in sequence_iterator
    merged = merge_fastq_reads(sequence_pair...)
    println(merged)
end

```


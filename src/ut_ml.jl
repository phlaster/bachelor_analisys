using FastaIO

function onehot_encode_dna(seq::AbstractString; revcomp=false)
    N = length(seq)
    # Инициализируем матрицу размером (4, N)
    encoding = zeros(Float32, 4, N)

    if revcomp
        # Словарь для комплементарной цепи
        comp_mapping = Dict{Char,UInt8}('A'=>0x1, 'C'=>0x2, 'G'=>0x3, 'T'=>0x4)
        # Перебираем последовательность в обратном порядке
        for (i, ch) in enumerate(Iterators.reverse(uppercase(seq)))
            if haskey(comp_mapping, ch)
                row = comp_mapping[ch]
                encoding[row, i] = 1.f0
            end
        end
    else
        # Словарь для прямой цепи
        mapping = Dict{Char,UInt8}('T'=>0x1, 'G'=>0x2, 'C'=>0x3, 'A'=>0x4)
        # Перебираем последовательность в прямом порядке
        for (i, ch) in enumerate(uppercase(seq))
            if haskey(mapping, ch)
                row = mapping[ch]
                encoding[row, i] = 1.f0
            end
        end
    end
    return encoding
end


function CDS_borders_in_GFF(gff_data, chromosomes::Vector{<:AbstractString}; strand="+")
    @assert !isempty(gff_data) "Empty gff data"
    @assert !isempty(chromosomes) "Empty chromosomes list"
    @assert strand in ["+", "-"] "wrong strand identifier: $strand"
    
    cds_starts = Vector{Float32}[]
    cds_stops = Vector{Float32}[]
    for chrom in chromosomes
        region_filter = filter_gff_region(; sequence_header=chrom, regiontype="region")

        chrom_length = gff_data |> region_filter |> ranges_from_GFF_records |> only |> maximum

        cds_filter = filter_gff_region(;
            sequence_header=chrom,
            regiontype="CDS",
            strand=strand
        )
        intervals = ranges_from_GFF_records(cds_filter(gff_data))
        

        chrom_starts = zeros(Float32, chrom_length)
        chrom_stops = zeros(Float32, chrom_length)

        if isempty(intervals)
            push!(cds_starts, chrom_starts)
            push!(cds_stops, chrom_stops)
            continue
        end

        starts = replace(
            getfield.(intervals, :start) .% chrom_length,
            0 => chrom_length
        )
        stops = replace(
            getfield.(intervals, :stop) .% chrom_length,
            0 => chrom_length
        )
        @assert length(starts) == length(stops) "start-stop count mismatch"
        chrom_starts[starts] .= 1.f0
        chrom_stops[stops] .= 1.f0
        
        push!(cds_starts, chrom_starts)
        push!(cds_stops, chrom_stops)
    end
    
    return strand == "+" ?
        (cds_starts, cds_stops) :
        (reverse!.(cds_starts), reverse!.(cds_stops))
end

function digitize_genome(genome_dir::T) where T <: AbstractString
    @assert isdir(genome_dir) "No such genomic directory"
    file_prefix = last(splitpath(genome_dir))
    fastaname = joinpath(genome_dir, "$file_prefix.fna")
    gffname = joinpath(genome_dir, "$file_prefix.gff")
    @assert all(isfile.([fastaname, gffname])) "Genomic files not found"
    
    fasta_data = readfasta(fastaname)
    fasta_strings = getindex.(fasta_data, 2)
    
    dna_matrix_pos = onehot_encode_dna.(fasta_strings; revcomp=false)
    dna_matrix_neg = onehot_encode_dna.(fasta_strings; revcomp=true)
    
    chromosomes = getindex.(fasta_data, 1) .|> split .|> first
    gff_data = open_gff(gffname)
    
    starts_pos, ends_pos = CDS_borders_in_GFF(gff_data, chromosomes; strand="+")
    starts_neg, ends_neg = CDS_borders_in_GFF(gff_data, chromosomes; strand="-")

    @assert all(==(length(fasta_strings)), [
        length(dna_matrix_pos),
        length(dna_matrix_neg),
        length(starts_pos),
        length(ends_pos),
        length(starts_neg),
        length(ends_neg)
    ]) "Dimentions do not match"

    return dna_matrix_pos, dna_matrix_neg, starts_pos, ends_pos, starts_neg, ends_neg
end
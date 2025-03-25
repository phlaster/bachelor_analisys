using FastaIO

function onehot_encode_dna(seq::AbstractString; revcomp::Bool=false)
    N = length(seq)
    encoding = zeros(Float32, 4, N)
    N == 0 && return encoding

    @inbounds if revcomp
        comp_mapping = Dict('T'=>1, 'G'=>2, 'C'=>3, 'A'=>4)
        for (i, ch) in enumerate(Iterators.reverse(uppercase(seq)))
            if haskey(comp_mapping, ch)
                row = comp_mapping[ch]
                encoding[row, i] = 1.f0
            end
        end
    else
        mapping = Dict('A'=>1, 'C'=>2, 'G'=>3, 'T'=>4)
        for (i, ch) in enumerate(uppercase(seq))
            if haskey(mapping, ch)
                row = mapping[ch]
                encoding[row, i] = 1.f0
            end
        end
    end
    return encoding
end

function onehot_decode_dna(encoding::Matrix{T}; revcomp::Bool=false) where T <: Number
    if size(encoding, 1) != 4
        error("Incorrect matrix dimentions: 4 rows expected, получено $(size(encoding, 1)).")
    end

    N = size(encoding, 2)
    decoded = Vector{Char}(undef, N)

    mapping = revcomp ? Dict(
        T[1, 0, 0, 0]=>'T',
        T[0, 1, 0, 0]=>'G',
        T[0, 0, 1, 0]=>'C',
        T[0, 0, 0, 1]=>'A'
    ) : Dict(
        T[1, 0, 0, 0]=>'A',
        T[0, 1, 0, 0]=>'C',
        T[0, 0, 1, 0]=>'G',
        T[0, 0, 0, 1]=>'T'
    )

    @inbounds @simd for j in 1:N
        col = encoding[:, j]
        if haskey(mapping, col)
            decoded[j] = mapping[col]
        else
            decoded[j] = 'N'
        end
    end

    return join(revcomp ? reverse!(decoded) : decoded)
end



function CDS_borders_in_GFF(gff_data, chrom_names, chrom_lengths; strand="+")
    @assert !isempty(gff_data) "Empty gff data"
    @assert !isempty(chrom_names) "Empty chromosomes list"
    @assert strand in ["+", "-"] "wrong strand identifier: $strand"
    
    cds_starts = Vector{Float32}[]
    cds_stops = Vector{Float32}[]
    for (chrom_name, chrom_length) in zip(chrom_names, chrom_lengths)
        # region_filter = filter_gff_region(; sequence_header=chrom_name, regiontype="region")
        # chrom_length = gff_data |> region_filter |> ranges_from_GFF_records |> only |> maximum

        cds_filter = filter_gff_region(;
            sequence_header=chrom_name,
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

function CDS_borders_one_side(gff_data, chrom_names, chrom_lengths, side::Symbol; strand="+")
    @assert !isempty(gff_data) "Empty gff data"
    @assert !isempty(chrom_names) "Empty chromosomes list"
    @assert strand in ["+", "-"] "wrong strand identifier: $strand, must be either \"+\" or \"-\""
    @assert side in [:starts, :stops] "wrong side symbol: $side, must be either :starts or :stops"
    
    cds_borders = [zeros(Float32, chrom_length) for chrom_length in chrom_lengths]
    for (i, chrom_name, chrom_length) in zip(1:999999, chrom_names, chrom_lengths)
        cds_filter = filter_gff_region(;
            sequence_header=chrom_name,
            regiontype="CDS",
            strand=strand
        )
        filtered_regions = cds_filter(gff_data)

        @inbounds if strand=="+"
            borders = side == :starts ? ranges_from_GFF_records(filtered_regions, :left) :
                              ranges_from_GFF_records(filtered_regions, :right)
            isempty(borders) && continue
            cycle_borders = replace(borders .% chrom_length, 0 => chrom_length)
            cds_borders[i][cycle_borders] .= 1.f0
        else
            borders = side == :starts ? ranges_from_GFF_records(filtered_regions, :right) :
                              ranges_from_GFF_records(filtered_regions, :left)
            isempty(borders) && continue
            cycle_borders = replace(borders .% chrom_length, 0 => chrom_length)
            cds_borders[i][cycle_borders] .= 1.f0
            reverse!(cds_borders[i])
        end
    end
    
    return cds_borders
end


function digitize_genome(genome_dir::T) where T <: AbstractString
    @assert isdir(genome_dir) "No such genomic directory"
    file_prefix = last(splitpath(genome_dir))
    fastaname = joinpath(genome_dir, "$file_prefix.fna")
    gffname = joinpath(genome_dir, "$file_prefix.gff")
    @assert all(isfile.([fastaname, gffname])) "Genomic files not found"
    
    fasta_data = readfasta(fastaname)
    fasta_strings = getindex.(fasta_data, 2)
    chrom_lengths = length.(fasta_strings)
    
    dna_matrix_pos = onehot_encode_dna.(fasta_strings; revcomp=false)
    dna_matrix_neg = onehot_encode_dna.(fasta_strings; revcomp=true)
    
    chrom_names = getindex.(fasta_data, 1) .|> split .|> first
    gff_data = open_gff(gffname)
    
    starts_pos, ends_pos = CDS_borders_in_GFF(gff_data, chrom_names, chrom_lengths; strand="+")
    starts_neg, ends_neg = CDS_borders_in_GFF(gff_data, chrom_names, chrom_lengths; strand="-")

    @assert all(==(length(fasta_strings)), [
        length(dna_matrix_pos),
        length(dna_matrix_neg),
        length(starts_pos),
        length(ends_pos),
        length(starts_neg),
        length(ends_neg)
    ]) "Dimentions do not match"
    
    res = (
        dna_pos=dna_matrix_pos,
        dna_neg=dna_matrix_neg,
        starts_pos=starts_pos,
        stops_pos=ends_pos,
        starts_neg=starts_neg,
        stops_neg=ends_neg
    )
    return res
end

function digitize_genome_one_side(genome_dir::T, side::Symbol=:starts) where T <: AbstractString
    @assert side in [:starts, :stops] "wrong side symbol: $side, must be either :starts or :stops"
    @assert isdir(genome_dir) "No such genomic directory"
    file_prefix = last(splitpath(genome_dir))
    fastaname = joinpath(genome_dir, "$file_prefix.fna")
    gffname = joinpath(genome_dir, "$file_prefix.gff")
    @assert all(isfile.([fastaname, gffname])) "Genomic files not found"
    
    fasta_data = readfasta(fastaname)
    fasta_strings = getindex.(fasta_data, 2)
    chrom_lengths = length.(fasta_strings)
    
    dna_pos = onehot_encode_dna.(fasta_strings; revcomp=false)
    dna_neg = onehot_encode_dna.(fasta_strings; revcomp=true)
    
    chrom_names = getindex.(fasta_data, 1) .|> split .|> first
    gff_data = open_gff(gffname)
    
    borders_pos = CDS_borders_one_side(gff_data, chrom_names, chrom_lengths, side; strand="+")
    borders_neg = CDS_borders_one_side(gff_data, chrom_names, chrom_lengths, side; strand="-")
    
    res = (
        dna_pos=dna_pos,
        dna_neg=dna_neg,
        borders_pos=borders_pos,
        borders_neg=borders_neg
    )
    return res
end
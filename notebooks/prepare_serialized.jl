### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 7f3c544c-0766-11f0-3da6-b3c5c6b653e3
begin
	using Pkg
	Pkg.activate(dirname(@__DIR__))
end

# ╔═╡ 06ea3164-06e5-4f0f-8dc1-7b70ee3e45ca
using FastaIO, SimplePlutoInclude

# ╔═╡ d2cef288-9178-47b8-945b-eca56ebe5cd2
using Serialization, ProgressLogging

# ╔═╡ ada341a4-9f03-476d-a70e-4268448c7e58
@plutoinclude "../src/utils.jl"

# ╔═╡ 246df942-64e6-4ce0-81f3-be42a2e6369f
accs = readdir("../DATA/genomes/genomes")

# ╔═╡ ab1a738a-206d-409a-b377-aa42dda21a7a
function split_intervals_by_remainder(vec_of_ranges::Vector{UnitRange{Int}})
    rems0 = UnitRange{Int64}[]
    rems1 = UnitRange{Int64}[]
    rems2 = UnitRange{Int64}[]

    for range in vec_of_ranges
        first_elem = first(range)
        rem = mod(first_elem, 3)
        if rem == 0
            push!(rems0, range)
        elseif rem == 1
            push!(rems1, range)
        else
            push!(rems2, range)
        end
    end

    return rems0, rems1, rems2
end

# ╔═╡ 5ba62567-4842-4763-a03a-0216e2c5d015
function anyintersect(ranges::Vector{UnitRange{Int}})
    n = length(ranges)
    for i in 1:(n-1)
        for j in (i+1):n
            low = max(first(ranges[i]), first(ranges[j]))
            high = min(last(ranges[i]), last(ranges[j]))
            if low <= high
				@warn "Intersection between $i-th ($(ranges[i])) and $j-th ($(ranges[j])) intervals"
                return true
            end
        end
    end
    return false
end

# ╔═╡ 8a8f5462-eff8-4d52-bd30-73b846bddd97
function all_ranges_to_feature_tensor!(tensor, all_ranges, n_max, strand_coord)
    num_sets = length(all_ranges)

	# позиции, стренды, фреймы, фичи
    @assert size(tensor) == (n_max, 2, num_sets, 3)
    
    for (set_idx, ranges) in enumerate(all_ranges)
        for interval in ranges
            a = first(interval)
            b = last(interval)
            
            if 1 <= a <= n_max
                tensor[a, strand_coord, set_idx, 1] = 0x1
            end
            if 1 <= b <= n_max
                tensor[b, strand_coord, set_idx, 3] = 0x1
            end
            for pos in (a+1):(b-1)
                if 1 <= pos <= n_max
                    tensor[pos, strand_coord, set_idx, 2] = 0x1
                end
            end
        end
    end
    return tensor
end


# ╔═╡ 6d90398e-5544-4c89-b899-2f66ec6912e3
function onehot_encode_dna(seq::AbstractString; revcomp=false)
    N = length(seq)
    encoding = zeros(UInt8, N, 4)

	if revcomp
		comp_mapping = Dict{Char,UInt8}('A'=>0x1, 'C'=>0x2, 'G'=>0x3, 'T'=>0x4)
	    for (i, ch) in enumerate(Iterators.reverse(uppercase(seq)))
	        if haskey(comp_mapping, ch)
	            col = comp_mapping[ch]
	            encoding[i, col] = 0x1
	        end
	    end
	else
		mapping = Dict{Char,UInt8}('T'=>0x1, 'G'=>0x2, 'C'=>0x3, 'A'=>0x4)
	    for (i, ch) in enumerate(uppercase(seq))
	        if haskey(mapping, ch)
	            col = mapping[ch]
	            encoding[i, col] = 0x1
	        end
	    end
	end
    return encoding
end

# ╔═╡ 389bddda-5acf-4f01-9da7-2bf717675f80
function CDS_borders_in_GFF(gff_data, chromosomes::Vector{<:AbstractString}; strand="+")
	@assert !isempty(gff_data) "Empty gff data"
	@assert !isempty(chromosomes) "Empty chromosomes list"
	@assert strand in ["+", "-"] "wrong strand identifier: $strand"
	
	cds_starts = Vector{UInt8}[]
	cds_stops = Vector{UInt8}[]
	for chrom in chromosomes
		region_filter = filter_gff_region(; sequence_header=chrom, regiontype="region")

		chrom_length = gff_data |> region_filter |> ranges_from_GFF_records |> only |> maximum

		cds_filter = filter_gff_region(;
			sequence_header=chrom,
			regiontype="CDS",
			strand=strand
		)
		intervals = ranges_from_GFF_records(cds_filter(gff_data))
		

		chrom_starts = zeros(UInt8, chrom_length)
		chrom_stops = zeros(UInt8, chrom_length)

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
		chrom_starts[starts] .= 0x1
		chrom_stops[stops] .= 0x1
		
		push!(cds_starts, chrom_starts)
		push!(cds_stops, chrom_stops)
	end
	
	return strand == "+" ?
		(cds_starts, cds_stops) :
		(reverse!.(cds_starts), reverse!.(cds_stops))
end

# ╔═╡ 678cdb2c-b59b-41a2-b5d4-08be245859f0
function digitize_genome(genome_dir)
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

# ╔═╡ 7dede4f3-0063-4578-9bf4-6e29474531c3
rand_dir = "../DATA/genomes/genomes/$(rand(accs))"

# ╔═╡ 3eaf6a63-182d-40f1-b96e-78324d403b4b
fixed_dir = "../DATA/genomes/genomes/$(accs[42])"

# ╔═╡ 7e432565-a843-4703-827b-9f09eea70c1b
# ╠═╡ disabled = true
#=╠═╡
@Threads.threads for (i, acc) in collect(enumerate(accs))
	dir = "../DATA/genomes/genomes/$acc"
	dg_name = "../DATA/genomes/serialized/$acc"
	isfile(dg_name) && continue
	try
		digitize_genome(dir)
		if mod(i, 100) == 0
			println(i)
		end
		# serialize(dg_name, dg)
	catch e
		@warn i, dir
		break
	end
end
  ╠═╡ =#

# ╔═╡ 085437f6-604d-40bc-bd34-7c7e024634e1
single =  digitize_genome("../DATA/genomes/genomes/GCF_000010825.1")

# ╔═╡ b2607609-9242-4609-9920-b334612cd191
single[3]

# ╔═╡ Cell order:
# ╠═7f3c544c-0766-11f0-3da6-b3c5c6b653e3
# ╠═06ea3164-06e5-4f0f-8dc1-7b70ee3e45ca
# ╠═ada341a4-9f03-476d-a70e-4268448c7e58
# ╠═246df942-64e6-4ce0-81f3-be42a2e6369f
# ╠═ab1a738a-206d-409a-b377-aa42dda21a7a
# ╠═5ba62567-4842-4763-a03a-0216e2c5d015
# ╠═8a8f5462-eff8-4d52-bd30-73b846bddd97
# ╠═6d90398e-5544-4c89-b899-2f66ec6912e3
# ╠═389bddda-5acf-4f01-9da7-2bf717675f80
# ╠═678cdb2c-b59b-41a2-b5d4-08be245859f0
# ╠═7dede4f3-0063-4578-9bf4-6e29474531c3
# ╠═3eaf6a63-182d-40f1-b96e-78324d403b4b
# ╠═d2cef288-9178-47b8-945b-eca56ebe5cd2
# ╠═7e432565-a843-4703-827b-9f09eea70c1b
# ╠═085437f6-604d-40bc-bd34-7c7e024634e1
# ╠═b2607609-9242-4609-9920-b334612cd191

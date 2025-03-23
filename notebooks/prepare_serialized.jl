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
using FastaIO, Serialization

# ╔═╡ 453350bb-1bc7-4b41-816c-8a9c69313271
include("../src/utils.jl")

# ╔═╡ 246df942-64e6-4ce0-81f3-be42a2e6369f
accs = readdir("../DATA/genomes/genomes")

# ╔═╡ ab1a738a-206d-409a-b377-aa42dda21a7a
function split_intervals_by_remainder(ranges::Vector{UnitRange{Int64}})
    rems0 = UnitRange{Int64}[]
    rems1 = UnitRange{Int64}[]
    rems2 = UnitRange{Int64}[]

    for range in ranges
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
function non_intersecting_intervals(ranges::Vector{UnitRange{Int64}})
    n = length(ranges)
    for i in 1:(n-1)
        for j in (i+1):n
            low = max(first(ranges[i]), first(ranges[j]))
            high = min(last(ranges[i]), last(ranges[j]))
            if low <= high
				println(ranges[i], "  ",ranges[j])
                return false
            end
        end
    end
    return true
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
function onehot_encode_dna(seq::AbstractString)
    N = length(seq)
    encoding = zeros(UInt8, N, 4)
    
    mapping = Dict{Char,UInt8}('A'=>0x1, 'G'=>0x2, 'T'=>0x3, 'C'=>0x4)
    
    for (i, ch) in enumerate(seq)
        if haskey(mapping, ch)
            col = mapping[ch]
            encoding[i, col] = 1
        end
    end
    
    return encoding
end

# ╔═╡ 678cdb2c-b59b-41a2-b5d4-08be245859f0
function assemble_tensor(genome_dir)
	file_prefix = last(splitpath(genome_dir))
	fastaname = joinpath(genome_dir, "$file_prefix.fna")
	gffname = joinpath(genome_dir, "$file_prefix.gff")
	
	fasta_data = readfasta(fastaname)
	gff_data = open_gff(gffname)

	full_dna_matrix = onehot_encode_dna(join(getindex.(fasta_data, 2)))

	chrom_names = getindex.(fasta_data, 1) .|> split .|> first
	@show length(chrom_names)
	tensors = []
	fasta_onehot = []
	for ch in chrom_names
		chrom_length = filter_gff_region(;
			sequence_header=string(ch), regiontype="region"
		)(gff_data) |> ranges_from_GFF_records |> only |> maximum
		@show chrom_length
		
		chrom_tenz = zeros(UInt8, chrom_length, 2, 3, 3)
		for (strand_dim, strand) in enumerate(["+", "-"])
			extraction_filter = filter_gff_region(;
				sequence_header=string(ch),
				regiontype="CDS",
				strand=strand
			)
			xtrcted_intrvls = ranges_from_GFF_records(extraction_filter(gff_data))
			
			frames_0_1_2 = xtrcted_intrvls |>
				(x->filter(y->mod(length(y),3)==0, x)) |>
				split_intervals_by_remainder

			all_ranges_to_feature_tensor!(
				chrom_tenz, frames_0_1_2, chrom_length, strand_dim
			)		
		end
		push!(tensors, chrom_tenz)
	end
	full_acc_tensor = cat(tensors..., dims=1)
	@show size(full_acc_tensor), size(full_dna_matrix)
	return full_acc_tensor, full_dna_matrix
end

# ╔═╡ 7dede4f3-0063-4578-9bf4-6e29474531c3
rand_dir = "../DATA/genomes/genomes/$(rand(accs))"

# ╔═╡ f2e3441a-f112-4c8a-8e71-de46504a09e2
@time random_genome_encoded = assemble_tensor(rand_dir);

# ╔═╡ Cell order:
# ╠═7f3c544c-0766-11f0-3da6-b3c5c6b653e3
# ╠═06ea3164-06e5-4f0f-8dc1-7b70ee3e45ca
# ╠═453350bb-1bc7-4b41-816c-8a9c69313271
# ╠═246df942-64e6-4ce0-81f3-be42a2e6369f
# ╠═ab1a738a-206d-409a-b377-aa42dda21a7a
# ╠═5ba62567-4842-4763-a03a-0216e2c5d015
# ╠═8a8f5462-eff8-4d52-bd30-73b846bddd97
# ╠═6d90398e-5544-4c89-b899-2f66ec6912e3
# ╠═678cdb2c-b59b-41a2-b5d4-08be245859f0
# ╠═7dede4f3-0063-4578-9bf4-6e29474531c3
# ╠═f2e3441a-f112-4c8a-8e71-de46504a09e2

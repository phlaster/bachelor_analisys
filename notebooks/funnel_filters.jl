### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ af402476-f84f-11ef-3dbf-d7dae3493197
using DataFrames, CSV, Chain

# ╔═╡ e8eb0a90-e8bf-4cdf-8d59-d8ac7954522e
metadata_file = "../DATA/tree_GTDB/bac120_metadata_r220.tsv"

# ╔═╡ 444dee8f-f21d-4717-b8be-9a2f8545ef3a
tree_table_file = "../DATA/tree_GTDB/clustered_tree.tsv"

# ╔═╡ 98b04859-92b9-483c-9b81-0006d214c2d0
g = (
	completeness_min=95,
	contamination_max=5,
	contig_count_max=200,
	trna_count_min=30,
	n50_min=50_000,
	sh_max=10,
	ambiguous_bases_max=0
)

# ╔═╡ 56249a14-b75d-4653-af1c-8328162a36fc
metadata_df = CSV.read(metadata_file, DataFrame);

# ╔═╡ e3b8d15f-3636-49ed-bdc6-61f34b4f8382
cluster_df = DataFrame(CSV.File(tree_table_file, delim='\t', header=true), ["accession", "cluster"]);

# ╔═╡ 8167863e-5b95-4239-b833-972da766a682
begin
	funnel = Pair{Int64, String}[]
@chain metadata_df begin
	@aside push!(funnel, nrow(_)=>"All")

	innerjoin(_, cluster_df, on=:accession)
	@aside push!(funnel, nrow(_)=>"In tree")
	
	# filter(r->r.gtdb_type_designation_ncbi_taxa_sources != "none", _)
	# @aside push!(funnel, nrow(_)=>"ncbi_taxa")
	
	# filter(r->occursin.(strip.(first.(split.(r.ncbi_organism_name)), Ref(['[', ']'])), r.gtdb_taxonomy), _ )
	# @aside push!(funnel, nrow(_)=>"ncbi=gtdb taxonomy")
	
	# filter(r->r.mimag_high_quality == "t", _)
	# @aside push!(funnel, nrow(_)=>"mimag hq")
	
	filter(r->r.ncbi_assembly_level in ["Complete Genome", "Chromosome"], _)
	@aside push!(funnel, nrow(_)=>"Genome+Chromosome")
	
	filter(r->r.ncbi_rrna_count != "none", _)
	@aside push!(funnel, nrow(_)=>">0 rRNA")
	
	filter(r->min(r.checkm2_completeness, r.checkm_completeness) ≥ g.completeness_min, _)
	@aside push!(funnel, nrow(_)=>"completeness>$(g.completeness_min)")

	filter(r->max(r.checkm2_contamination, r.checkm_contamination) ≤ g.contamination_max, _)
	@aside push!(funnel, nrow(_)=>"contamination<$(g.contamination_max)")
	
	# filter(r->r.contig_count ≤ g.contig_count_max, _)
	# @aside push!(funnel, nrow(_)=>"contigs<$(g.contig_count_max)")
	
	filter(r->r.trna_count ≥ g.trna_count_min, _)
	@aside push!(funnel, nrow(_)=>"tRNA>$(g.trna_count_min)")
	
	filter(r->r.n50_contigs ≥ g.n50_min, _)
	@aside push!(funnel, nrow(_)=>"n50>$(g.n50_min)")
	
	# filter(r->r.checkm_strain_heterogeneity ≤ g.sh_max, _)
	# @aside push!(funnel, nrow(_)=>"SH<$(g.sh_max)")
	
	filter(r->r.ambiguous_bases == g.ambiguous_bases_max, _)
	@aside push!(funnel, nrow(_)=>"amb bases ≤ $(g.ambiguous_bases_max)")
end
	funnel
end

# ╔═╡ ce71bae7-7533-463e-abce-12e96660dcad
print(funnel)

# ╔═╡ c03a4c9a-5d96-4c89-b7dc-9053fddb094d
metadata_df.gtdb_type_designation_ncbi_taxa_sources |> Set

# ╔═╡ e5ae5cdc-ef03-4709-b092-965b77d89666
function generate_sankey_code(pairs::Vector{Pair{Int64, String}})
    if isempty(pairs)
        error("Input vector is empty!")
    end
    
    lines = String[]
    first_pair = pairs[1]
    # push!(lines, "Input [$(first_pair.first)] $(first_pair.second)")
    for i in 1:(length(pairs) - 1)
        current = pairs[i]
        nxt = pairs[i+1]
        pass_value = nxt.first
        drop_value = current.first - nxt.first        
        push!(lines, "$(current.second) [$(pass_value)] $(nxt.second)")
        if drop_value > 0
            push!(lines, "$(current.second) [$(drop_value)] -$(round(100*drop_value/(pass_value+drop_value); digits=2))%")
        end
    end
    last_pair = pairs[end]
    push!(lines, "$(last_pair.second) [$(last_pair.first)] Output")
    
    return join(lines, "\n")
end

# ╔═╡ 1fdd8271-d8f5-4df5-9262-6c7c7cb786c8
generate_sankey_code(funnel) |> print

# ╔═╡ 992df2b7-8761-45bd-8b1c-846acf9712d2
let
	funnel = Pair{Int64, String}[]
@chain metadata_df begin
	@aside push!(funnel, nrow(_)=>"All")

	innerjoin(_, cluster_df, on=:accession)
	@aside push!(funnel, nrow(_)=>"In tree")
	
	# filter(r->r.gtdb_type_designation_ncbi_taxa_sources != "none", _)
	# @aside push!(funnel, nrow(_)=>"ncbi_taxa")
	
	# filter(r->occursin.(strip.(first.(split.(r.ncbi_organism_name)), Ref(['[', ']'])), r.gtdb_taxonomy), _ )
	# @aside push!(funnel, nrow(_)=>"ncbi=gtdb taxonomy")
	
	filter(r->r.mimag_high_quality == "t", _)
	@aside push!(funnel, nrow(_)=>"mimag hq")
	
	# filter(r->r.ncbi_assembly_level in ["Complete Genome", "Chromosome"], _)
	# @aside push!(funnel, nrow(_)=>"Genome+Chromosome")
	
	# filter(r->r.ncbi_rrna_count != "none", _)
	# @aside push!(funnel, nrow(_)=>">0 rRNA")
	
	# filter(r->min(r.checkm2_completeness, r.checkm_completeness) ≥ g.completeness_min, _)
	# @aside push!(funnel, nrow(_)=>"completeness>$(g.completeness_min)")

	# filter(r->max(r.checkm2_contamination, r.checkm_contamination) ≤ g.contamination_max, _)
	# @aside push!(funnel, nrow(_)=>"contamination<$(g.contamination_max)")
	
	# filter(r->r.contig_count ≤ g.contig_count_max, _)
	# @aside push!(funnel, nrow(_)=>"contigs<$(g.contig_count_max)")
	
	# filter(r->r.trna_count ≥ g.trna_count_min, _)
	# @aside push!(funnel, nrow(_)=>"tRNA>$(g.trna_count_min)")
	
	# filter(r->r.n50_contigs ≥ g.n50_min, _)
	# @aside push!(funnel, nrow(_)=>"n50>$(g.n50_min)")
	
	# filter(r->r.checkm_strain_heterogeneity ≤ g.sh_max, _)
	# @aside push!(funnel, nrow(_)=>"SH<$(g.sh_max)")
	
	# filter(r->r.ambiguous_bases == g.ambiguous_bases_max, _)
	# @aside push!(funnel, nrow(_)=>"amb bases ≤ $(g.ambiguous_bases_max)")
end
	funnel
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
Chain = "8be319e6-bccf-4806-a6f7-6fae938471bc"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"

[compat]
CSV = "~0.10.15"
Chain = "~0.6.0"
DataFrames = "~1.7.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.3"
manifest_format = "2.0"
project_hash = "c61938a0cc6596fadbb68c5d89183b0b526b45b5"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

[[deps.Chain]]
git-tree-sha1 = "9ae9be75ad8ad9d26395bf625dea9beac6d519f1"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.6.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fb61b4812c49343d7ef0b533ba982c46021938a6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.7.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "2ec417fc319faa2d768621085cc1feebbdee686b"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.23"

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

    [deps.FilePathsBase.weakdeps]
    Mmap = "a63ad114-7e13-5084-954f-fe012c677804"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.InlineStrings]]
git-tree-sha1 = "6a9fde685a7ac1eb3495f8e812c5a7c3711c2d5e"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.3"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "725421ae8e530ec29bcbdddbe91ff8053421d023"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╠═af402476-f84f-11ef-3dbf-d7dae3493197
# ╠═e8eb0a90-e8bf-4cdf-8d59-d8ac7954522e
# ╠═444dee8f-f21d-4717-b8be-9a2f8545ef3a
# ╠═98b04859-92b9-483c-9b81-0006d214c2d0
# ╠═56249a14-b75d-4653-af1c-8328162a36fc
# ╠═e3b8d15f-3636-49ed-bdc6-61f34b4f8382
# ╠═8167863e-5b95-4239-b833-972da766a682
# ╠═ce71bae7-7533-463e-abce-12e96660dcad
# ╠═c03a4c9a-5d96-4c89-b7dc-9053fddb094d
# ╠═e5ae5cdc-ef03-4709-b092-965b77d89666
# ╠═1fdd8271-d8f5-4df5-9262-6c7c7cb786c8
# ╠═992df2b7-8761-45bd-8b1c-846acf9712d2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

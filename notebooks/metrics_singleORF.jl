### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ cbe0c970-d345-11ef-2576-b5ddb0e0a7b0
using GFF3, SimplePlutoInclude, DataFrames, Format, PrettyTables

# ╔═╡ 8daa248e-0608-49b7-aac7-25b40927529d
Base.show(io::IO, x::Float64) = print(io, format( x, precision=3 ))

# ╔═╡ 7cbf2a02-bc90-4355-be87-5949643fffa4
md"""
```jl
println(io, " seqid: ", hasseqid(record)       ? seqid(record) : "<...>")
println(io, "source: ", hassource(record)      ? source(record) : "<...>")
println(io, "  type: ", hasfeaturetype(record) ? featuretype(record) : "<...>")
println(io, " start: ", hasseqstart(record)    ? seqstart(record) : "<...>")
println(io, "   end: ", hasseqend(record)      ? seqend(record) : "<...>")
println(io, " score: ", hasscore(record)       ? score(record) : "<...>")
println(io, "strand: ", hasstrand(record)      ? strand(record) : "<...>")
println(io, " phase: ", hasphase(record)       ? phase(record) : "<...>")
```
"""

# ╔═╡ bfb5c23b-69b8-4ed4-bfb6-0df00da3b29b
@plutoinclude "../src/utils.jl"

# ╔═╡ 94f0d899-b0e9-4cfe-907c-23744c91abfe
reference_file = "../DATA/GCF_004306555.1_ASM430655v1_genomic.gff";

# ╔═╡ cd51d4b0-63be-4491-9f36-196a68d37a7d
helixer_dir = "../DATA/helixer_results/"

# ╔═╡ b6eabca1-2c13-48f4-9e86-3714b32b7910
helixer_files = filter(endswith(".gff"), readdir(helixer_dir, join=true))

# ╔═╡ 356685d0-380d-4b3c-a678-ea296f10434b
chroms = Set([GFF3.seqid(x) for x in open_gff(reference_file)]) |> collect |> sort

# ╔═╡ 9298083d-6bb4-4a55-8403-742ae00c568d
chrom_1 = chroms[1]

# ╔═╡ 313a2490-6fc3-40d1-a3c7-c647aa02bdde
confmatrix(ref) = x -> ConfusionMTR(ref, x);

# ╔═╡ dbc1f3bb-8933-439b-b7ff-de72dd46c578
pipeline = ranges_from_GFF_records ∘ filter_gff_region(chrom_1; regiontype="CDS", strand="+", phase=0) ∘ open_gff;

# ╔═╡ 484011a1-b4e0-4cf8-a0e6-56df96fd96d7
reference_ranges = pipeline(reference_file);

# ╔═╡ 87454cc7-3605-4dd0-9d15-ce02ead175dc
helixir_ranges = pipeline.(helixer_files);

# ╔═╡ 7e9735d5-e627-45c4-854f-8b355d68240b
ranges_all = merge_ranges(helixir_ranges);

# ╔═╡ 0dde3c6a-dc59-4c25-b0b0-1c2ad160b7b0
CMs = helixir_ranges .|> confmatrix(reference_ranges)

# ╔═╡ d2a6048c-9e77-45bb-9fdd-de3c9b2398a3
CM_comb = confmatrix(reference_ranges)(ranges_all);

# ╔═╡ bcc63fd8-12d2-458d-b37e-c3477dd48779
let
	data = []
	for (name, cm) in zip(helixer_files, CMs)
	    push!(data, [
			first(split(basename(name), '.')), cm.prec, cm.rec, cm.f1, cm.fdr
		])
	end
	push!(data, ["combined", CM_comb.prec, CM_comb.rec, CM_comb.f1, CM_comb.fdr])
	
	pretty_table(permutedims(hcat(data...)); header=["Sample", "Precision", "Recall", "F1", "FDR"])
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Format = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
GFF3 = "af1dc308-cb6b-11e8-32f0-31192efa90f6"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
SimplePlutoInclude = "6f00a2c5-ea4a-46bf-9183-91b7b57a087f"

[compat]
DataFrames = "~1.7.0"
Format = "~1.3.7"
GFF3 = "~0.2.3"
PrettyTables = "~2.4.0"
SimplePlutoInclude = "~0.2.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.3"
manifest_format = "2.0"
project_hash = "c6e4119d635129cf6231d68aa6c63f7ea519c2c2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Automa]]
deps = ["TranscodingStreams"]
git-tree-sha1 = "ef9997b3d5547c48b41c7bd8899e812a917b409d"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.4"

[[deps.BGZFStreams]]
deps = ["CodecZlib"]
git-tree-sha1 = "3aca54d25f8c30056577aa37ea68184da68df685"
uuid = "28d598bf-9b8f-59f1-b38c-5a06b4a0f5e6"
version = "0.3.2"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BioGenerics]]
deps = ["TranscodingStreams"]
git-tree-sha1 = "017562e86afcd2a6a2a9220606a40b54604887c9"
uuid = "47718e42-2ac5-11e9-14af-e5595289c2ea"
version = "0.1.5"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

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

[[deps.FASTX]]
deps = ["Automa", "BioGenerics", "PrecompileTools", "ScanByte", "StringViews", "TranscodingStreams"]
git-tree-sha1 = "5c4b85ab007cb55d38fc8249ddfab6bf2f48cf06"
uuid = "c2308a5c-f048-11e8-3e8a-31650f418d12"
version = "2.1.2"

    [deps.FASTX.extensions]
    BioSequencesExt = "BioSequences"

    [deps.FASTX.weakdeps]
    BioSequences = "7e6ae17a-c86d-528c-b3b9-7f778a29fe59"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GFF3]]
deps = ["Automa", "BGZFStreams", "BioGenerics", "FASTX", "GenomicFeatures", "Indexes", "TranscodingStreams", "URIParser"]
git-tree-sha1 = "ace846d50b5e1c0ea2cbb13032dc11165ebf8383"
uuid = "af1dc308-cb6b-11e8-32f0-31192efa90f6"
version = "0.2.3"

[[deps.GenomicFeatures]]
deps = ["BioGenerics", "DataStructures", "IntervalTrees"]
git-tree-sha1 = "51b0906aab4a9ae1a4095c27994361a582c5ebe1"
uuid = "899a7d2d-5c61-547b-bef9-6698a8d05446"
version = "2.1.0"

[[deps.Indexes]]
deps = ["BGZFStreams", "BioGenerics", "GenomicFeatures", "TranscodingStreams"]
git-tree-sha1 = "275bce824b40fd2e70358e0a652ba1b34172f240"
uuid = "4ffb77ac-cb80-11e8-1b35-4b78cc642f6d"
version = "0.1.3"

[[deps.InlineStrings]]
git-tree-sha1 = "45521d31238e87ee9f9732561bfee12d4eebd52d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.2"

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

[[deps.IntervalTrees]]
git-tree-sha1 = "dc3b97bb5c9cb7c437f74027309f2c2f09a82aaf"
uuid = "524e6230-43b7-53ae-be76-1e9e4d08d11b"
version = "1.1.0"

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

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
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

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OrderedCollections]]
git-tree-sha1 = "12f1439c4f986bb868acda6ea33ebc78e19b95ad"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.7.0"

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

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "fea870727142270bdf7624ad675901a1ee3b4c87"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.7.1"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "d49e35f413186528f1d7cc675e67d0ed16fd7800"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.4.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SimplePlutoInclude]]
deps = ["Dates"]
git-tree-sha1 = "a98dee6bdc63f648d92b28348dcb9b59c2881fbf"
uuid = "6f00a2c5-ea4a-46bf-9183-91b7b57a087f"
version = "0.2.0"

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
git-tree-sha1 = "a6b1675a536c5ad1a60e5a5153e1fee12eb146e3"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.0"

[[deps.StringViews]]
git-tree-sha1 = "ec4bf39f7d25db401bcab2f11d2929798c0578e5"
uuid = "354b36f9-a18e-4713-926e-db85100087ba"
version = "1.3.4"

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

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.URIParser]]
deps = ["Unicode"]
git-tree-sha1 = "53a9f49546b8d2dd2e688d216421d050c9a31d0d"
uuid = "30578b45-9adc-5946-b283-645ec420af67"
version = "0.4.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

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
# ╠═cbe0c970-d345-11ef-2576-b5ddb0e0a7b0
# ╠═8daa248e-0608-49b7-aac7-25b40927529d
# ╟─7cbf2a02-bc90-4355-be87-5949643fffa4
# ╠═bfb5c23b-69b8-4ed4-bfb6-0df00da3b29b
# ╠═94f0d899-b0e9-4cfe-907c-23744c91abfe
# ╠═cd51d4b0-63be-4491-9f36-196a68d37a7d
# ╠═b6eabca1-2c13-48f4-9e86-3714b32b7910
# ╠═356685d0-380d-4b3c-a678-ea296f10434b
# ╠═9298083d-6bb4-4a55-8403-742ae00c568d
# ╠═313a2490-6fc3-40d1-a3c7-c647aa02bdde
# ╠═dbc1f3bb-8933-439b-b7ff-de72dd46c578
# ╠═484011a1-b4e0-4cf8-a0e6-56df96fd96d7
# ╠═87454cc7-3605-4dd0-9d15-ce02ead175dc
# ╠═7e9735d5-e627-45c4-854f-8b355d68240b
# ╠═0dde3c6a-dc59-4c25-b0b0-1c2ad160b7b0
# ╠═d2a6048c-9e77-45bb-9fdd-de3c9b2398a3
# ╠═bcc63fd8-12d2-458d-b37e-c3477dd48779
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

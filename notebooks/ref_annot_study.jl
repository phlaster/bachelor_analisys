### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 0625cb2c-d831-11ef-14f7-8b07feee0d27
using GFF3, SimplePlutoInclude, Format, DataFrames, PrettyTables

# ╔═╡ 295e2401-1658-4d1f-af11-7f4763a59ca3
Base.show(io::IO, x::Float64) = print(io, format( x, precision=3 ))

# ╔═╡ 6cd26414-76be-4197-814d-8db8fde46a6a
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

# ╔═╡ e2ed4dad-b1ac-4f06-bd3f-b91d5d78fbe3
@plutoinclude "../src/utils.jl"

# ╔═╡ 1bf3d2a3-6034-491b-971b-5a5981fea30b
begin
	reference_file = "../DATA/GCF_004306555.1_ASM430655v1_genomic.gff";
	fungi_file = "../DATA/helixer_results/fungi_helixer.gff"
	invertebrate_file = "../DATA/helixer_results/invertebrate_helixer.gff"
	land_plant_file = "../DATA/helixer_results/land_plant_helixer.gff"
	vertebrate_file = "../DATA/helixer_results/vertebrate_helixer.gff"
end

# ╔═╡ db6bda7d-a82b-4431-b912-e16b9fd3d86d
chrom_1 = "NZ_SILH01000001.1";

# ╔═╡ 8e5bb5f1-6f2f-47a8-85ee-e2c14c6c3c78
md"""
# Найденные регионы
| Тип региона          | Описание                                                                                            |
|----------------------|:----------------------------------------------------------------------------------------------------|
| **gene**             | Основной функциональный элемент, кодирующий РНК или белок.                                          |
| **CDS**              | Кодирующая последовательность (Coding Sequence) — часть гена, кодирующая белок.                     |
| **exon**             | Экзон — кодирующая часть гена (актуально для эукариот).                                             |
| **tRNA**             | Транспортная РНК (transfer RNA) — участвует в трансляции.                                           |
| **pseudogene**       | Псевдоген — нефункциональный остаток гена.                                                          |
| **riboswitch**       | Рибопереключатель — участок РНК, регулирующий экспрессию генов.                                     |
| **rRNA**             | Рибосомная РНК (ribosomal RNA) — компонент рибосом.                                                 |
| **region**           | Общий регион (например, хромосома или плазмида).                                                    |
| **tmRNA**            | Трансфер-мессенджерная РНК (transfer-messenger RNA) — участвует в разрешении проблем с трансляцией. |
| **ncRNA**            | Некодирующая РНК (non-coding RNA) — функциональная РНК, не кодирующая белок.                        |
| **sequence_feature** | Общий термин для любого функционального элемента последовательности.                                |
| **SRP_RNA**          | РНК сигнального пептида (Signal Recognition Particle RNA) — участвует в транспорте белков.          |
| **RNase\_P\_RNA**      | РНК компонент фермента RNase P — участвует в процессинге тРНК.                                      |
"""

# ╔═╡ d0648b01-718f-4456-9172-641ed64e7392
md"""
## Исходная аннотация
"""

# ╔═╡ 0f7433e4-2026-4982-a474-a9673f770d13
feature_strand_table(open_gff(reference_file), chrom_1)

# ╔═╡ a917583d-406d-4e89-a7aa-f8f73b2eabf0
rand_locus(open_gff(reference_file), chrom_1, "exon"; strand="+")

# ╔═╡ 0a2c74dc-c2dd-49f8-8a89-2ea8907af6fa
md"""
## fungi
"""

# ╔═╡ 6c8d2287-7cb2-4df0-a1c7-9dd95c201dae
feature_strand_table(open_gff(fungi_file), chrom_1)

# ╔═╡ 679361e2-1638-4ccd-8829-70ce7c35e9e4
md"""
## land_plant
"""

# ╔═╡ e81a3d04-8bc6-4a3d-93c8-4e0a41601376
feature_strand_table(open_gff(land_plant_file), chrom_1)

# ╔═╡ aaf0e7d7-151e-4390-8f17-30ec6c6202ea
md"""
## vertebrate
"""

# ╔═╡ 695d0d63-fb00-4472-9bce-252e3b93a016
feature_strand_table(open_gff(vertebrate_file), chrom_1)

# ╔═╡ 137b831f-b863-4472-9a56-98e97011da55
md"""
## invertebrate
"""

# ╔═╡ ba08e658-58f9-40a3-ac71-e734b34b3e33
feature_strand_table(open_gff(invertebrate_file), chrom_1)

# ╔═╡ c992af98-3dfa-4504-b216-d6d1a3b156b6
rand_locus(open_gff(invertebrate_file), chrom_1, "exon")

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
# ╠═0625cb2c-d831-11ef-14f7-8b07feee0d27
# ╠═295e2401-1658-4d1f-af11-7f4763a59ca3
# ╟─6cd26414-76be-4197-814d-8db8fde46a6a
# ╠═e2ed4dad-b1ac-4f06-bd3f-b91d5d78fbe3
# ╠═1bf3d2a3-6034-491b-971b-5a5981fea30b
# ╠═db6bda7d-a82b-4431-b912-e16b9fd3d86d
# ╟─8e5bb5f1-6f2f-47a8-85ee-e2c14c6c3c78
# ╟─d0648b01-718f-4456-9172-641ed64e7392
# ╠═0f7433e4-2026-4982-a474-a9673f770d13
# ╠═a917583d-406d-4e89-a7aa-f8f73b2eabf0
# ╟─0a2c74dc-c2dd-49f8-8a89-2ea8907af6fa
# ╠═6c8d2287-7cb2-4df0-a1c7-9dd95c201dae
# ╟─679361e2-1638-4ccd-8829-70ce7c35e9e4
# ╠═e81a3d04-8bc6-4a3d-93c8-4e0a41601376
# ╟─aaf0e7d7-151e-4390-8f17-30ec6c6202ea
# ╠═695d0d63-fb00-4472-9bce-252e3b93a016
# ╟─137b831f-b863-4472-9a56-98e97011da55
# ╠═ba08e658-58f9-40a3-ac71-e734b34b3e33
# ╠═c992af98-3dfa-4504-b216-d6d1a3b156b6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

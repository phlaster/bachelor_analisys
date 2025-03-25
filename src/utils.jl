const UT_TECHNICAL = joinpath(@__DIR__, "ut_technical.jl")
const UT_CONFUSION = joinpath(@__DIR__, "ut_confusion.jl")
const UT_GFF_RANGES = joinpath(@__DIR__, "ut_gff_ranges.jl")
const UT_ML = joinpath(@__DIR__, "ut_ml.jl")

include(UT_TECHNICAL)
include(UT_CONFUSION)
include(UT_GFF_RANGES)
include(UT_ML)






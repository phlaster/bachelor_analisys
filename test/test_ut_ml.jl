#!/usr/bin/env julia
const PROJECT_DIR = dirname(@__DIR__)
const UTILS_FILE = joinpath(PROJECT_DIR, "src", "utils.jl")
const PLOTS_FILE = joinpath(PROJECT_DIR, "src", "plots.jl")

using Pkg
Pkg.activate(PROJECT_DIR, io=devnull)
include(UTILS_FILE)
using Test

@testset "true_distances tests" begin

    # Test 1: Example given in the problem.
    v = [false, false, true, true, false, false, true, false, true]
    expected = [1, 3, 2]
    @test true_distances(v) == expected

    # Test 2: Single true value returns an empty array.
    v = [false, false, true, false]
    @test true_distances(v) == []

    # Test 3: No true values returns an empty array.
    v = [false, false, false, false]
    @test true_distances(v) == []

    # Test 4: Multiple trues separated by one false.
    v = [true, false, true, false, true]
    # Indices are [1, 3, 5] so distances should be [2, 2]
    @test true_distances(v) == [2, 2]

    # Test 5: All true values.
    v = [true, true, true, true]
    # Indices: [1,2,3,4] -> distances [1,1,1]
    @test true_distances(v) == [1, 1, 1]

    # Test 6: Trues at the beginning and end.
    v = [true, false, false, false, true]
    # Indices: [1,5] -> distance [4]
    @test true_distances(v) == [4]

    # Test 7: Two trues adjacent.
    v = [true, true]
    @test true_distances(v) == [1]

end

@testset "false_positive_stats tests" begin
    # Test 1: Simple case with one FP to match
    ground_truth = [false, false, true, true, false, false, true, false, true]
    predicted    = [false, false, false, true, false, true, false, false, true]
    # Explanation:
    # - Real positives in ground_truth: positions [3, 4, 7, 9]
    # - Predicted true: positions [4, 6, 9]
    #   => positions 4 and 9 are direct matches.
    #   => Only FP is at position 6. The only unmatched real positive is position 7.
    #   => Signed distance = 7 - 6 = 1.
    @test false_positive_stats(ground_truth, predicted) == [1.0]

    # Test 2: Multiple FP’s with nontrivial distances.
    # Ground truth has real positives at positions: 2, 6, 10, 14
    # Predicted has direct matches at positions: 6, 14 and false positives at positions: 3, 11, 15
    # After matching direct ones, unmatched real: [2, 10]
    # FP candidates: 3, 11, 15
    # Compute candidate pairs:
    # For FP 3: 
    #   distance to 2: 2-3 = -1 (abs 1); to 10: 7 (abs 7)
    # For FP 11:
    #   distance to 2: 2-11 = -9 (abs 9); to 10: 10-11 = -1 (abs 1)
    # For FP 15:
    #   distance to 2: -13 (abs 13); to 10: -5 (abs 5)
    # Sorted by abs distance: (1, FP3, real2, -1), (1, FP11, real10, -1), (5, FP15, real10, -5), ...
    # First, FP3 is paired with real2 giving -1, and FP11 is paired with real10 giving -1.
    # FP15 is left without a candidate (its candidate real10 is already matched).
    ground_truth = falses(15)
    predicted    = falses(15)
    for i in [2,6,10,14]
      ground_truth[i] = true
    end
    for i in [3,6,11,14,15]
      predicted[i] = true
    end
    # Direct matches at positions 6 and 14; unmatched real positives: 2 and 10; FP candidates: 3, 11, 15.
    # Expected result: [-1.0, -1.0] (order not important if algorithm collects by closeness)
    result = false_positive_stats(ground_truth, predicted)
    sort!(result)
    @test result == [-1.0, -1.0]

    # Test 3: No false positives (all predicted positives are exact matches)
    ground_truth = [true, false, true, false]
    predicted    = [true, false, true, false]
    @test false_positive_stats(ground_truth, predicted) == Float64[]

    # Test 4: All predictions false, so no FP.
    ground_truth = [true, true, false]
    predicted    = [false, false, false]
    @test false_positive_stats(ground_truth, predicted) == Float64[]

    # Test 5: More FP candidates than unmatched real positives.
    # Real positives: positions 4 and 10.
    # Predicted: additional FPs: positions 1, 3, 8, 12.
    # Direct match: assume position 4 is a direct match, so unmatched real: [10].
    # Then we have FP candidates: 1, 3, 8, 12 (removing 4 if it were predicted true, but we set only one direct match).
    # Let’s construct:
    ground_truth = [false, false, false, true, false, false, false, false, false, true, false, false]
    predicted    = [true, false, true, true, false, false, false, true, false, false, false, true]
    # Here, direct match: position 4.
    # FP candidates: 1, 3, 8, 12.
    # Unmatched real: [10]
    # We expect pairing FP with real 10:
    #  For FP 8: distance = 10-8 = 2 (abs=2)
    #  For FP 12: distance = 10-12 = -2 (abs=2)
    #  For FP 1: distance = 9 (abs=9)
    #  For FP 3: distance = 7 (abs=7)
    # The algorithm will choose one candidate with minimum absolute distance. Note that there is a tie between FP8 and FP12.
    # Because we sorted candidates in order of increasing abs (and then by order in the candidate array), the one encountered first will be taken.
    # The exact candidate pairing may depend on the candidate order. We need to check that exactly one pairing is made and its distance is either 2 or -2.
    result = false_positive_stats(ground_truth, predicted)
    @test length(result) == 1
    @test (result[1] == 2.0 || result[1] == -2.0)

    # Test 6: All predicted positives are false positives.
    # Real positives: positions [5, 9]. None of these are predicted.
    # Predicted: some other indices are true.
    ground_truth = [false, false, false, false, true, false, false, false, true, false]
    predicted    = [true, false, true, false, false, false, true, false, false, true]
    # So direct matches: none.
    # FP candidates: positions [1,3,7,10].
    # Unmatched real: positions [5, 9].
    # Candidate distances:
    # FP1: [5-1=4; 9-1=8] -> (4,8)
    # FP3: [5-3=2; 9-3=6] -> (2,6)
    # FP7: [5-7=-2; 9-7=2] -> (2,2)
    # FP10: [5-10=-5; 9-10=-1] -> (5, -1)
    # Sorted by absolute distances:
    # FP3 to real5: abs(2)=2, FP7 to real5: abs(2)=2, FP7 to real9: abs(2)=2, FP10 to real9: abs(1)=1, ...
    # Actually, check FP10 to real9: 9-10 = -1, abs=1 is the smallest.
    # So the algorithm should match FP10 with real9 giving -1.
    # Then the next pair among remaining candidates will be chosen from FP3, FP7 with unmatched real left [5].
    # Among those FP3 to real5: 5-3 = 2 (abs 2) and FP7 to real5: -2 (abs 2). Tie-breaker order may choose FP3 pairing.
    # So expected distances: [-1.0, 2.0] or [-1.0, -2.0] depending on candidate ordering.
    result = false_positive_stats(ground_truth, predicted)
    @test length(result) == 2
    # Check that result contains -1.0 and one of 2.0 or -2.0.
    @test -1.0 in result
    @test (2.0 in result || -2.0 in result)

    # Test 7: Vectors of different length should throw an assertion error.
    ground_truth = [true, false, true]
    predicted = [true, false]
    @test_throws AssertionError false_positive_stats(ground_truth, predicted)
end
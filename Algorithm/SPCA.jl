# SPCA optimizer
# Arguments:
# 1. Input CSV. One header row expected, no row header (all rows have N entries)
# 2. Output CSV. Feature loadings are appended to an open file.
# 3. Parameter K. The number of variables in the support of each PC.
# 4. Parameter D. The iteration count for the multiComponents loop.
# 5. Parameter searchCap. Halt optimization early because we generally have:
# searchCap << N nCr K.

include("utilities.jl")
include("branchAndBound.jl")
include("multiComponents.jl")

using CSV, DataFrames, LinearAlgebra, Random, Tables

LinearAlgebra.BLAS.set_num_threads(2)

EIGEN_GAP_ARG = 6
RANDOM_SEED_ARG = 7

if length(ARGS) >= RANDOM_SEED_ARG
    Random.seed!(parse(UInt64, ARGS[RANDOM_SEED_ARG]))
end

K = parse(Int64, ARGS[3])
D = parse(Int64, ARGS[4])
searchCap = parse(Int64, ARGS[5])

eigenGap = parse(Float64, ARGS[EIGEN_GAP_ARG])

csvFile = CSV.File(ARGS[1])
genes = csvFile.names
Sigma = csvFile |> Tables.matrix

# print((;zip(genes, repeat([]))))
# print((; zip(genes, repeat([Float64[]], length(genes)))...))
# empty_output = DataFrame(
#     (; zip(genes, repeat([Float64[]], length(genes)))...)
# )
empty_output = DataFrame([name => [] for name in genes])
outputFile = ARGS[2]
CSV.write(outputFile, empty_output; quotestrings=true)

# Run multiOptimalSPCA
myprob = problem(zeros(size(Sigma)), Sigma)
for i = 1:D
    # data = repeat([0.0], length(genes))
    ~, xVal = branchAndBound(
        myprob, K, outputFlag=0, gap=eigenGap, searchCap=searchCap,
        id=ARGS[RANDOM_SEED_ARG]
    )
    myprob.Sigma = projectMatrix(myprob.Sigma, xVal)
    CSV.write(
        outputFile,
        # Tables.table(data'),
        Tables.table(xVal'),
        append=true
    )
end

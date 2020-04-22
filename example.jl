include("src/index.jl")
using .GameTree
using Distributions
using Plots


results = []
    plot_errors = []
    plot_samples = []
    runs = 10
    ϵ = 0.1
    δ = 0.1
    for i in 1:runs
        tree = createTree(minnode, 3, 10)
        @time mu = FindTopWinner2!(ϵ, 0.1, tree)
        error = abs(mu[tree][1] - tree.val)
        samples = sum(map(x -> x[2], values(mu)))
        push!(plot_errors, error)
        push!(plot_samples, samples)
        push!(results, samples)
    end
    histogram(results)

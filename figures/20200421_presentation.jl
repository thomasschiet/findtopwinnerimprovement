include("src/index.jl")
using .GameTree
using Distributions
using Plots

# 1. Solving c(δ)
K = 4
    δ = 0.5
    f(x) = exp(-x)*(2*exp(1)*x/K)^(K/2)
    x_start = -0.8
    x_end = 6
    xs = x_start:0.01:x_end
    ys = f.(xs)
    p = plot(xs, ys, xlims=(x_start, x_end), label="\$\\delta(c) = e^{-c}(2ec/K)^{K/2}\$")
    plot!(x -> δ, label = "\$\\delta  = 0.5\$")
    xlabel!("c")
    title!("\$K=4\$")
    # title!(string("\$\\delta\\textrm{ as function of c. } K = ", K, "\$"))
    plot!(size=(600/4*3, 400/4*3), dpi=200)
    savefig(p, "fig.png")
    p

# 2. Solving c(δ)
K = 4
    δ = 0.5
    f(x) = exp(-x)*(2*exp(1)*x/K)^(K/2)
    x_start = 2
    x_end = 6
    xs = x_start:0.01:x_end
    ys = f.(xs)
    p = plot(xs, ys, xlims=(x_start, x_end), label="\$\\delta(c) = e^{-c}(2ec/K)^{K/2}\$")
    plot!(x -> δ, label = "\$\\delta  = 0.5\$")
    xlabel!("c")
    title!("\$K=4\$")
    # title!(string("\$\\delta\\textrm{ as function of c. } K = ", K, "\$"))
    plot!(size=(600/4*3, 400/4*3), dpi=200)
    savefig(p, "fig2.png")
    p

# 3.1 Results FTW1
results = []
    plot_errors = []
    plot_samples = []
    runs = 250
    ϵ = 0.1
    δ = 0.1
    depth = 3
    width = 10
    type = minnode
    for i in 1:runs
        tree = createTree(type, depth, width)
        # treecopy = deepcopy(tree)
        @time mu = FindTopWinner!(ϵ, δ, tree)
        error = abs(mu[tree][1] - tree.val)
        samples = sum(map(x -> x[2], values(mu)))
        push!(plot_errors, error)
        push!(plot_samples, samples)
        push!(results, (samples, error))
    end
    println("Mean error:", mean(plot_errors))
    println("# outside error bound: ", 1 - count(x -> x < ϵ, plot_errors)/length(plot_errors))
    println("Mean samples: ", mean(plot_samples))
    scatter(plot_samples, plot_errors, label = "FindTopWinner")

# 3.2 Results FTW2
results_2 = []
    plot_errors_2 = []
    plot_samples_2 = []
    global tree
    for i in 1:runs
        tree = createTree(type, depth, width)
        @time mu = FindTopWinner2!(ϵ, δ, tree)
        error = abs(mu[tree][1] - tree.val)
        samples = sum(map(x -> x[2], values(mu)))
        push!(plot_errors_2, error)
        push!(plot_samples_2, samples)
        push!(results_2, (samples, error))
    end
    println("Mean error:", mean(plot_errors_2))
    println("# outside error bound: ", 1 - count(x -> x < ϵ, plot_errors_2)/length(plot_errors_2))
    println("Mean samples: ", mean(plot_samples_2))
    scatter!(plot_samples_2, plot_errors_2, label = "FindTopWinner improvement")

xlims!((0, 2*1e6))
xaxis!("# oracle calls")
yaxis!("Absolute error")
plot!(size=(600/4*3, 400/4*3), dpi=200)
savefig("results.png")
# FindTopWinner
# Mean error:0.0091738248385439
# #outside error bound: 0.0
# Mean samples: 1.27401658e6

# FindTopWinner2
# Mean error:0.008677946831329107
# # outside error bound: 0.0
# Mean samples: 719649.14



# FindTopWinner prunes more in first round, but less in later rounds

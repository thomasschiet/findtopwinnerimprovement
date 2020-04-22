MuDict = Dict{GameTreeNode, Tuple{Float64, Int}}

function getMu(tree)::MuDict
    getMu(tree, MuDict())
end

function getMu(tree, dict::MuDict)::MuDict
    dict[tree] = (0, 0)
    for child in tree.children
        getMu(child, dict)
    end
    dict
end

function FindTopWinner!(eps, delta, tree::GameTreeNode)
    root = tree
    eps_m = 1
    delta_m = delta/count_leaves(tree)

    mu = getMu(tree)

    rounds = ceil(log2(MathConstants.e/eps))
    for m in 1:rounds
        eps_m /= 2
        delta_m /= 2
        if length(tree.children) <= 1
            break
        end
        mu[root] = EstimateValues!(root, eps_m, delta_m, mu)
        Prune!(root, eps_m, mu)
        # println("Epoch: ", m, " #leaves: ", count_leaves(root))
    end

    mu
end

function EstimateValues!(node::GameTreeNode, eps_m::Float64, delta_m::Float64, mu::MuDict)::Tuple{Float64, Int}
    if node.type == leaf
        totalCalls = Int(ceil((1/(2*eps_m^2)) * log(2/delta_m)))
        newCalls = totalCalls - mu[node][2]
        newMean = mean(rand(node.oracle, newCalls)) * newCalls + mu[node][1] * mu[node][2]
        newMean /= totalCalls
        return (newMean, totalCalls)
    end

    for child in node.children
        mu[child] = EstimateValues!(child, eps_m, delta_m, mu)
    end

    if node.type == maxnode
        return maximum(map(child -> mu[child], node.children))
    else
        return minimum(map(child -> mu[child], node.children))
    end
end

function Prune!(node, eps_m, mu)
    prunedChildren = filter(child -> abs(mu[node][1] - mu[child][1]) > 2eps_m, node.children)
    filter!(child -> child ∉ prunedChildren, node.children)

    for child in node.children
        Prune!(child, eps_m, mu)
    end
end

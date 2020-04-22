using Distributions

@enum NodeType minnode=1 maxnode=2 leaf=3

struct GameTreeNode
    val::Union{Float64, Missing}
    type::NodeType
    children::Vector{GameTreeNode}
    oracle::Union{Missing, Sampleable}
end

GameTreeNode(val, type, children) = GameTreeNode(val, type, children, missing)

"""
Count number of nodes in a tree `tree`
"""
function count_nodes(tree::GameTreeNode)::Int
    1 + reduce(+, map(count_nodes, tree.children), init = 0)
end

"""
Count number of leaves in a tree `tree`
"""
function count_leaves(tree::GameTreeNode)::Int
    if tree.type == leaf
        1
    else
        sum(map(count_leaves, tree.children))
    end
end

"""
Creates a game tree of heigth `heigth` and outdegree `width`.
The oracles have a normal distribution with mean uniformly in [0, 1] and variance 1.
"""
function createTree(rootType::NodeType, height::Int, width::Int)::GameTreeNode
    @assert height >= 0
    @assert width > 0
    @assert !(height == 0 && rootType != leaf)

    childrenType = if height == 1
        leaf
    elseif rootType == maxnode
        minnode
    else
        maxnode
    end
    children = width * (rootType != leaf)

    children = map(x -> createTree(childrenType, height - 1, width), zeros(children))

    if rootType == leaf
        val = rand()
        oracle = Normal(val, 1)
    elseif rootType == maxnode
        val = maximum(map(x -> x.val, children))
        oracle = missing
    else
        val = minimum(map(x -> x.val, children))
        oracle = missing
    end

    GameTreeNode(val, rootType, children, oracle)
end

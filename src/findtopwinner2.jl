using Random, Distributions
using Plots
using Test

function FindTopWinner2!(eps, delta, tree::GameTreeNode)
    root = tree
    eps_m = 1
    delta_m = delta/count_leaves(tree)

    mu = getMu(tree)

    for m in 1:(ceil(log2(MathConstants.e/eps)))
        eps_m /= 2
        delta_m /= 2
        if length(tree.children) <= 1
            break
        end
        mu[root] = EstimateValues2!(root, eps_m, delta_m, mu)
        Prune2!(root, eps_m, delta_m, mu)
        # println("Epoch: ", m, " #leaves: ", count_leaves(root))
    end

    mu
end

function EstimateValues2!(node::GameTreeNode, eps_m::Float64, delta_m::Float64, mu::MuDict)::Tuple{Float64, Int}
    if node.type == leaf
        totalCalls = Int(ceil((1/(2*eps_m^2)) * log(2/delta_m)))
        newCalls = totalCalls - mu[node][2]
        newMean = mean(rand(node.oracle, newCalls)) * newCalls + mu[node][1] * mu[node][2]
        newMean /= totalCalls
        return (newMean, totalCalls)
    end

    for child in node.children
        mu[child] = EstimateValues2!(child, eps_m, delta_m, mu)
    end

    if node.type == maxnode
        return maximum(map(child -> mu[child], node.children))
    else
        return minimum(map(child -> mu[child], node.children))
    end
end

function Prune2!(node, eps_m, delta_m, mu)
    if node.type == leaf
        return
    end
    @assert node.type != leaf
    @assert length(node.children) > 0

    if length(node.children) > 1
        prunedChildren = filter(child -> shouldPruneChild(mu, node, child, eps_m, delta_m), node.children)

        if length(prunedChildren) >= length(node.children)
            @assert length(prunedChildren) < length(node.children)
        end

        filter!(child -> child ∉ prunedChildren, node.children)
    end


    for child in node.children
        Prune2!(child, eps_m, delta_m, mu)
    end
end


function shouldPruneChild(mu::MuDict, parent::GameTreeNode, child::GameTreeNode, eps_m::Float64, delta_m)::Bool
    other_children = filter(x -> x != child, parent.children)
    if length(other_children) == 0
        return false
    end
    calls = ceil((1/(2*eps_m^2)) * log(2/delta_m))
    rho_squared = 2*getC(count_leaves(parent), delta_m) / calls
    @assert rho_squared > 0

    if parent.type == maxnode
        virtual_node = GameTreeNode(missing, maxnode, other_children, missing)

        # Make UCB_c(αρ^2) by restricting the UCB to [0, ρ^2] and then scale it down to [0, 1]
        UCB_c = scale(restrictto(UCB_on_tree(mu, child), rho_squared), 1/rho_squared)
        # Make LCB_u((1-α)ρ^2) by restricting the LCB to [0, ρ^2] and reflecting it about 1/2 ρ^2 and then scale it down to [0, 1]
        LCB_u_reflected = scale(restrict_and_reflect(LCB_on_tree(mu, virtual_node), rho_squared), 1/rho_squared)

        !(maximumvalue(UCB_c - LCB_u_reflected, (0., 1.)) > 0)
    else
        virtual_node = GameTreeNode(missing, minnode, other_children, missing)
        LCB_c = scale(restrictto(LCB_on_tree(mu, child), rho_squared), 1/rho_squared)
        UCB_u_reflected = scale(restrict_and_reflect(UCB_on_tree(mu, virtual_node), rho_squared), 1/rho_squared)

        !(maximumvalue(UCB_u_reflected - LCB_c, (0., 1.)) > 0)
    end
end

# Calculate LCB for a [type]node, given the LCB's of its children
function LCB(type::NodeType, LCB_cs::Vector{PiecewiseLinear})
    @assert length(LCB_cs) > 0
    if type == maxnode
        LCB_max_node(LCB_cs)
    elseif type == minnode
        LCB_min_node(LCB_cs)
    end
end

# Calculate LCB for a maxnode, given the LCB's of its children
function LCB_max_node(LCB_cs)

    # Collect all pieces of the LCB's in one array
    # Add index of each LCB to the pieces
    i = 1
    pieces = [];
    for LCB_c in LCB_cs
        append!(pieces, map(piece -> (piece.rho, piece.val, piece.slope, i), LCB_c))
        i += 1
    end
    sort!(pieces, by = piece -> piece[2], rev = true)

    current_rho = 0
    current_cost_per_child = zeros(length(LCB_cs))
    @assert length(pieces) > 0
    result = Piece[];
    prev_val = pieces[1][2]
    for piece in pieces
        leap = prev_val - piece[2]
        prev_val = piece[2]
        current_cost = sum(current_cost_per_child)
        current_rho += leap * sum(current_cost_per_child)
        current_cost_per_child[piece[4]] += -1.0/piece[3]
        next_cost = sum(current_cost_per_child)
        push!(result, Piece(current_rho, piece[2], -1/next_cost))
    end

    # Sanity check
    for p in result
        @assert p.rho >= 0
        @assert p.slope <= 0
    end
    result
end

function LCB_min_node(LCB_cs::Vector{PiecewiseLinear})
    @assert length(LCB_cs) > 0

    if length(LCB_cs) > 2
        y = LCB_min_node(LCB_cs[1:2])
        x = [[y]; LCB_cs[3:end]];
        return LCB_min_node(x)
    end

    if length(LCB_cs) == 1
        return LCB_cs[1]
    end

    pieces = [
        sort(LCB_cs[1], by = piece -> piece.rho),
        sort(LCB_cs[2], by = piece -> piece.rho),
    ]

    minimumfunction(pieces[1], pieces[2])
end

function UCB(type::NodeType, UCB_cs::Vector{PiecewiseLinear})
    if type == minnode
        UCB_min_node(UCB_cs)
    elseif type == maxnode
        return UCB_max_node(UCB_cs)
    else
        error("Got $type")
    end
end

function UCB_min_node(UCB_cs)
    i = 1
    pieces = [];
    for UCB_c in UCB_cs
        append!(pieces, map(piece -> (piece.rho, piece.val, piece.slope, i), UCB_c))
        i += 1
    end
    sort!(pieces, by = piece -> piece[2])
    current_rho = 0
    current_cost_per_child = zeros(length(UCB_cs))

    result = Piece[];
    prev_val = pieces[1][2]
    for piece in pieces
        leap = piece[2] - prev_val
        prev_val = piece[2]
        current_cost = sum(current_cost_per_child)
        current_rho += leap * sum(current_cost_per_child)
        current_cost_per_child[piece[4]] += 1.0/piece[3]
        next_cost = sum(current_cost_per_child)
        push!(result, Piece(current_rho, piece[2], 1/next_cost))
    end
    result
end

function UCB_max_node(UCB_cs::Vector{PiecewiseLinear})
    @assert length(UCB_cs) > 0

    if length(UCB_cs) > 2
        y = UCB_max_node(UCB_cs[1:2])
        x = [[y]; UCB_cs[3:end]];
        return UCB_max_node(x)
    end

    if length(UCB_cs) == 1
        return UCB_cs[1]
    end

    pieces = [
        sort(UCB_cs[1], by = piece -> piece.rho),
        sort(UCB_cs[2], by = piece -> piece.rho),
    ]

    maximumfunction(pieces[1], pieces[2])
end

function UCB_on_tree(mu::MuDict, node::GameTreeNode)::PiecewiseLinear
    if node.type == leaf
        [Piece(0.0, mu[node][1], 1.0)]
        # [UCBPiece(mu[node][1])]
    elseif node.type == maxnode
        UCB(maxnode, map(x -> UCB_on_tree(mu, x), node.children))
    elseif node.type == minnode
        UCB(minnode, map(x -> UCB_on_tree(mu, x), node.children))
    end
end

function LCB_on_tree(mu::MuDict, node::GameTreeNode)::PiecewiseLinear
    if node.type == leaf
        [Piece(0.0, mu[node][1], -1.0)]
    elseif node.type == maxnode
        LCB(maxnode, map(x -> LCB_on_tree(mu, x), node.children))
    elseif node.type == minnode
        LCB(minnode, map(x -> LCB_on_tree(mu, x), node.children))
    end
end

function boundsontree(mu::MuDict, node::GameTreeNode)
    UCB_on_tree(mu, node), LCB_on_tree(mu, node)
end

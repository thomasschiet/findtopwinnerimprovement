include("src/index.jl")
using .GameTree
using Distributions
using Test

function deterministicTest()
    deterministicTest(FindTopWinner!)
    deterministicTest(FindTopWinner2!)
end
function deterministicTest(alg)
    # Trivial tree, should always return value 1.0
    tree_template = GameTreeNode(1.0, maxnode, GameTreeNode[
        GameTreeNode(1.0, leaf, [], Bernoulli(1.0)),
        GameTreeNode(0.5, leaf, [], Bernoulli(0.5)),
        GameTreeNode(0.0, leaf, [], Bernoulli(0.0)),
    ], missing)
        tree = deepcopy(tree_template)
        mu = alg(0.1, 0.1, tree)
        @test mu[tree][1] == 1
        @test tree_template.val == 1
        @assert count_leaves(tree) == 1

    # Trivial tree, should always return value 0.0
    tree_template = GameTreeNode(0.0, minnode, GameTreeNode[
        GameTreeNode(1.0, leaf, [], Bernoulli(1.0)),
        GameTreeNode(0.5, leaf, [], Bernoulli(0.5)),
        GameTreeNode(0.0, leaf, [], Bernoulli(0.0)),
    ], missing)
        tree = deepcopy(tree_template)
        mu = alg(0.1, 0.1, tree)
        @test mu[tree][1] == 0
        @test tree_template.val == 0
        @assert count_leaves(tree) == 1

    # Trivial tree, should always return value 0.0
    tree_template = GameTreeNode(0.0, minnode, GameTreeNode[
        GameTreeNode(1.0, leaf, [], Bernoulli(1.0)),
        GameTreeNode(0.5, leaf, [], Bernoulli(0.5)),
        GameTreeNode(0.0, leaf, [], Bernoulli(0.0)),
        createTree(maxnode, 2, 10),
    ], missing)
        tree = deepcopy(tree_template)
        mu = alg(0.1, 0.1, tree)
        @test mu[tree][1] == 0
        @test tree_template.val == 0
end

deterministicTest()

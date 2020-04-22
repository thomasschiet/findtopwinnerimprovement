module GameTree

export Piece, PiecewiseLinear
export getPlotPiecewise, plotPiecewise, plotPiecewise!
export UCB, UCB_on_tree
export LCB, LCB_on_tree
export findIntersections
export boundsontree
export evalpiecewise
export restrictto
export restrict_and_reflect
export findIntersections
export scale
export maximumvalue, minimumvalue
export prune!

export NodeType, GameTreeNode, count_nodes, count_leaves, createTree
export minnode, maxnode, leaf
export FindTopWinner!
export FindTopWinner2!

include("invertfunction.jl")
include("piecewise.jl")
include("gametree.jl")

include("findtopwinner.jl")
include("findtopwinner2.jl")

end

import Base.+, Base.-

struct Piece
    rho::Float64
    val::Float64
    slope::Float64
end
PiecewiseLinear = Vector{Piece}
# UCBPiece(val::Float64) = Piece(0.0, val, 1.0)
# LCBPiece(val::Float64) = Piece(0.0, val, -1.0)

"""
Makes a `PiecewiseLinear` variable callable.

f = [Piece(0.0, 0.0, 1.0)]
f(2.0) # 2.0
"""
function (f::PiecewiseLinear)(rho::Float64)::Float64
    evalpiecewise(f, rho)
end

"""
Return f(ρ)
"""
function evalpiecewise(f::PiecewiseLinear, rho)::Float64
    if rho < 0
        throw(DomainError(rho, "ρ should be nonnegative"))
    end
    interval_end = f[1].rho
    for i in 1:length(f)
        interval_start = interval_end
        interval_end = i >= length(f) ? Inf : f[i+1].rho

        if interval_start <= rho < interval_end
            return f[i].val + f[i].slope * (rho - f[i].rho)
        end
    end
end

function +(x::PiecewiseLinear, y::PiecewiseLinear)::PiecewiseLinear
    both = sort([x; y], by = piece -> piece.rho)
    current_rho = 0
    current_slope = 0
    result = Piece[]
    for piece in both
        if piece.rho > current_rho
            if length(result) > 0
                previous_rho = result[end].rho
                previous_val = result[end].val
            else
                previous_rho = 0
                previous_val = x[1].val + y[1].val
            end
            current_val = previous_val + current_slope * (current_rho - previous_rho)
            push!(result, Piece(current_rho, current_val, current_slope))
            current_rho = piece.rho
        end
        current_slope += piece.slope
    end
    if length(result) > 0
        previous_rho = result[end].rho
        previous_val = result[end].val
    else
        previous_rho = 0
        previous_val = x[1].val + y[1].val
    end
    current_val = previous_val + current_slope * (current_rho - previous_rho)
    push!(result, Piece(current_rho, current_val, current_slope))

    sort(result, by = piece -> piece.rho)
end

-(y::PiecewiseLinear) = map(piece -> Piece(piece.rho, -piece.val, -piece.slope), y)
-(x::PiecewiseLinear, y::PiecewiseLinear) = x + -y

"""
Restrict `f` to [0, `rho`]
"""
function restrictto(f::PiecewiseLinear, rho::Float64)::PiecewiseLinear
    @assert length(f) > 0
    @assert f[1].rho == 0
    result = filter(piece -> piece.rho < rho, f)
    if result[end] != rho
        slope = result[end].slope

        val = result[end].val + (rho - result[end].rho) * slope
        push!(result, Piece(rho, val, slope))
    end

    sort(result, by = piece -> piece.rho)
end

"""
Scale `f` by `alpha`
"""
function scale(f::PiecewiseLinear, alpha::Float64)::PiecewiseLinear
    map(piece -> Piece(piece.rho * alpha, piece.val, piece.slope / alpha), f)

end

"""
Restrict `f` to [0, `rho`] and then return `x -> f(rho - x)`
"""
function restrict_and_reflect(f::PiecewiseLinear, rho::Float64)::PiecewiseLinear
    sort(map(piece -> Piece(rho - piece.rho, piece.val, -piece.slope), restrictto(f, rho)), by = piece -> piece.rho)
end

"""
Returns lists `xs`, `ys` of points where the slope of `f` changes.
"""
function getPlotPiecewise(f::PiecewiseLinear)
    sorted_points = sort(f, by = piece -> piece.rho)
    xs = map(piece -> piece.rho, sorted_points)
    ys = map(piece -> piece.val, sorted_points)

    inc = (maximum([xs[end] * 0.5, 2]))
    x_last = xs[end] + inc
    y_last = sorted_points[end].val + inc * sorted_points[end].slope

    append!(xs, x_last)
    append!(ys, y_last)

    xs, ys
end

function plotPiecewise(f::PiecewiseLinear; kwargs...)
    xs, ys = getPlotPiecewise(f)
    plot(xs, ys; kwargs...)
end

function plotPiecewise!(f::PiecewiseLinear; kwargs...)
    xs, ys = getPlotPiecewise(f)
    plot!(xs, ys; kwargs...)
end

"""
Returns the maximum value attained by `f` on the interval `[x, y]`.
I.e. max_{ρ ∈ [x, y]} f(ρ).
"""
function maximumvalue(f::PiecewiseLinear, (x, y)::Tuple{Float64, Float64})
    in_interval = filter(piece -> x < piece.rho < y, f)
    max_center_val = if length(in_interval) > 0
        maximum(piece -> piece.val, in_interval)
    else
        -Inf
    end
    max(max_center_val, f(x), f(y))
end

"""
Returns the minimum value attained by `f` on the interval `[x, y]`.
I.e. min_{ρ ∈ [x, y]} f(ρ).
"""
function minimumvalue(f::PiecewiseLinear, (x, y)::Tuple{Float64, Float64})
    in_interval = filter(piece -> x < piece.rho < y, f)
    min_center_val = if length(in_interval) > 0
        minimum(piece -> piece.val, in_interval)
    else
        Inf
    end
    min(min_center_val, f(x), f(y))
end

"""
Returns (x, y) of intersection of `a` and `b` if it exists. 
If intersection is not unique, it returns (a_rho, a_val).
"""
function intersection(a::Piece, b::Piece)::Union{Tuple{Float64, Float64}, Nothing}
    if a.slope == b.slope
        if (b.val - b.slope * b.rho) == (a.val - a.slope * a.rho)
            (a.rho, a.val)
        end
    else
        rho = (b.val - b.slope * b.rho) - (a.val - a.slope * a.rho)
        rho /= a.slope - b.slope

        val = a.slope * rho + (a.val - a.slope * a.rho)

        (rho, val)
    end
end


"""
Find all intersections of piecewise linear functions `main_line` and `other_line`
"""
function findIntersections(type, main_line, other_line)
    # find intersections
    intersections = Piece[]
    for i in 1:length(main_line)
        piece = main_line[i]
        next_piece_rho = if i == length(main_line)
            Inf
        else
            main_line[i + 1].rho
        end
        for j in 1:length(other_line)
            other_piece = other_line[j]
            next_other_piece_rho = if j == length(other_line)
                Inf
            else
                other_line[j + 1].rho
            end
            if other_piece.rho > next_piece_rho
                continue
            end

            point = intersection(piece, other_piece)
            if point == nothing || point[1] > next_piece_rho || point[1] < piece.rho || point[1] > next_other_piece_rho
                continue
            end

            # intersection exists and is in interval
            if type == :max
                slope = maximum([other_piece.slope, piece.slope])
            else
                slope = minimum([other_piece.slope, piece.slope])
            end
            intersection_piece = Piece(point[1], point[2], slope)
            push!(intersections, intersection_piece)

        end
    end
    intersections
end



struct PieceSegment
    rho_begin::Float64
    rho_end::Float64
    val::Float64
    slope::Float64
end

"""
Add start and end of interval of pieces. (a `Piece` only contains the start)
"""
function addSegments(a::PiecewiseLinear)::Vector{PieceSegment}
    sorted = sort(a, by = piece -> piece.rho)
    result = PieceSegment[]
    for i in 1:length(a)
        piece = a[i]
        next_rho = if i == length(a)
            Inf
        else
            a[i + 1].rho
        end
        push!(result, PieceSegment(piece.rho, next_rho, piece.val, piece.slope))
    end
    result
end

"""
Returns x -> min(a(x), b(x))
"""
function minimumfunction(a::PiecewiseLinear, b::PiecewiseLinear)::PiecewiseLinear
    intersections = addSegments(findIntersections(:min, a, b))
    a = addSegments(a)
    b = addSegments(b)

    rhos = map(piece -> piece.rho_begin, [a; b; intersections])
    sort!(rhos)
    unique!(rhos)
    result = Piece[]
    # rho = rhos[1]
    for rho in rhos
        p = filter(piece -> piece.rho_begin <= rho < piece.rho_end, [a; b])
        lowest_val = minimum(piece -> piece.val + piece.slope*(rho - piece.rho_begin), p)
        filter!(piece -> piece.val + piece.slope*(rho - piece.rho_begin) == lowest_val, p)
        lowest_slope = minimum(piece -> piece.slope, p)
        if length(result) == 0 || result[end].slope != lowest_slope
            push!(result, Piece(rho, lowest_val, lowest_slope))
        end
    end
    result
end

"""
Returns x -> max(a(x), b(x))
"""
maximumfunction(a::PiecewiseLinear, b::PiecewiseLinear)::PiecewiseLinear = -minimumfunction(-a, -b)

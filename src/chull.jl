"""
    convex_hull(ps::Vector{T}) where T<:Tuple{Number, Number}

Compute the convex hull of a finite set of 2D points using Andrew’s monotone-chain algorithm
in O(n log n) time. Given a vector of `(x, y)` tuples, returns the vertices of the convex hull
in counterclockwise order, starting from the left-most, bottom-most point.

References
- https://doi.org/10.1016/0020-0190(79)90072-3
- https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
"""

function convex_hull(ps::Vector{T}) where T<:Tuple{Number,Number}
    isempty(ps) && Memento.error(_LOGGER, "At least one point is necessary for evaluating the convex hull")
    pts = sort(ps)
    length(pts) ≤ 2 && return pts   # trivial hulls

    # cross‐product: >0 ⇒ left turn
    cross(a::T, b::T, c::T) = (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1])

    lower = T[]; upper = T[]
    for p in pts
        while length(lower) ≥ 2 && cross(lower[end-1], lower[end], p) ≤ 0
            pop!(lower)
        end
        push!(lower, p)
    end
    for p in reverse(pts)
        while length(upper) ≥ 2 && cross(upper[end-1], upper[end], p) ≤ 0
            pop!(upper)
        end
        push!(upper, p)
    end

    pop!(lower); pop!(upper)   # remove duplicate endpoints
    vcat(lower, upper)         # full CCW hull
end

# """
#     _orientation(x::Tuple{<:Number, <:Number}, y::Tuple{<:Number, <:Number}, z::Tuple{<:Number, <:Number})

# Utility function for `convex_hull`. Given an ordered triplet, this function returns if three
# points are collinear, oriented in clock-wise or anticlock-wise direction. 
# """
# function _orientation(
#     x::Tuple{<:Number,<:Number}, 
#     y::Tuple{<:Number,<:Number}, 
#     z::Tuple{<:Number,<:Number}
#     )
#     a = (y[2] - x[2]) * (z[1] - y[1]) - (y[1] - x[1]) * (z[2] - y[2])
#     QCO.is_zero(a) ? 0 : (a > 0 ? 1 : 2)
# end

# """
#     _lt_filter(a::Tuple{<:Number, <:Number}, b::Tuple{<:Number, <:Number})

# Utility function for sorting step in `convex_hull`. Given two points, `a` and `b`, this function 
# returns true if `a` has larger polar angle (counterclock-wise direction) than `b` w.r.t. first point `chull_p1`. 
# """
# function _lt_filter(
#     a::Tuple{<:Number, <:Number}, 
#     b::Tuple{<:Number, <:Number}
#     )
#     o = QCO._orientation(chull_p1, a, b)
#     o == 0 ? LA.norm(collect(chull_p1 .- b))^2 >= LA.norm(collect(chull_p1 .- a))^2 : o == 2
# end
# """
#     convex_hull(points::Vector{T}) where T<:Tuple{Number, Number}

# Graham's scan algorithm to compute the convex hull of a finite set of `n` points in a plane 
# with time complexity `O(n*log(n))`. Given a vector of tuples of co-ordinates, this function returns a 
# vector of tuples of co-ordinates which form the convex hull of the given set of points. 

# Sources: https://doi.org/10.1016/0020-0190(72)90045-2
#          https://en.wikipedia.org/wiki/Graham_scan 
# """
# function convex_hull(points::Vector{T}) where T<:Tuple{Number, Number}
#     num_points = size(points)[1]

#     if num_points == 0
#         Memento.error(_LOGGER, "Atleast one point is necessary for evaluating the convex hull")
#     elseif num_points <= 2 
#         return points
#     end

#     min_y = points[1][2]
#     min = 1
    
#     # Bottom-most or else the left-most point when tie
#     for i in 2:num_points
#         y = points[i][2]
#         if (y < min_y) || (y == min_y && points[i][1] < points[min][1])
#             min_y = y
#             min = i
#         end
#     end
    
#     # Placing the bottom/left-most point at first position
#     points[1], points[min] = points[min], points[1]

#     global chull_p1 = points[1]

#     points = Base.sort(points, lt = QCO._lt_filter)

#     # If two or more points are collinear with chull_p1, remove all except the farthest from chull_p1
#     arr_size = 2
#     for i in 2:num_points
#         while (i < num_points) && (QCO._orientation(chull_p1, points[i], points[i+1]) == 0)
#             i += 1
#         end
#         points[arr_size] = points[i]
#         arr_size += 1
#     end

#     if arr_size <= 3
#         return []
#     end

#     chull = []
#     append!(chull, [points[1]])
#     append!(chull, [points[2]])
#     append!(chull, [points[3]])
    
#     for i in 4:(arr_size-1)
#         while (size(chull)[1] > 1) && (QCO._orientation(chull[end - 1], chull[end], points[i]) != 2)
#             pop!(chull)
#         end
#         append!(chull, [points[i]])
#     end
    
#     return chull
# end
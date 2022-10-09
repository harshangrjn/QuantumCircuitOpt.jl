"""
    _orientation(x::Tuple{<:Number, <:Number}, y::Tuple{<:Number, <:Number}, z::Tuple{<:Number, <:Number})

Utility function for `convex_hull`. Given an ordered triplet, this function returns if three
points are collinear, oriented in clock-wise or anticlock-wise direction. 
"""
function _orientation(x::Tuple{<:Number, <:Number}, 
                      y::Tuple{<:Number, <:Number}, 
                      z::Tuple{<:Number, <:Number})
    a = ((y[2] - x[2]) * (z[1] - y[1]) - (y[1] - x[1]) * (z[2] - y[2]))

    if isapprox(a, 0, atol=1e-6)
        return 0 # collinear
    elseif a > 0 
        return 1 # clock-wise
    else 
        return 2 # counterclock-wise
    end
end

"""
    _lt_filter(a::Tuple{<:Number, <:Number}, b::Tuple{<:Number, <:Number})

Utility function for sorting step in `convex_hull`. Given two points, `a` and `b`, this function 
returns true if `a` has larger polar angle (counterclock-wise direction) than `b` w.r.t. first point `chull_p1`. 
"""
function _lt_filter(a::Tuple{<:Number, <:Number}, b::Tuple{<:Number, <:Number})
    o = QCO._orientation(chull_p1, a, b)

    if o == 0
        if LA.norm(collect(chull_p1 .- b))^2 >= LA.norm(collect(chull_p1 .- a))^2
            return true
        else
            return false
        end
    elseif o == 2
        return true
    else
        return false
    end
end

"""
    convex_hull(points::Vector{T}) where T<:Tuple{Number, Number}

Graham's scan algorithm to compute the convex hull of a finite set of `n` points in a plane 
with time complexity `O(n*log(n))`. Given a vector of tuples of co-ordinates, this function returns a 
vector of tuples of co-ordinates which form the convex hull of the given set of points. 

Sources: https://doi.org/10.1016/0020-0190(72)90045-2
         https://en.wikipedia.org/wiki/Graham_scan 
"""
function convex_hull(points::Vector{T}) where T<:Tuple{Number, Number}
    
    num_points = size(points)[1]

    if num_points == 0
        Memento.error(_LOGGER, "Atleast one point is necessary for evaluating the convex hull")
    elseif num_points <= 2 
        return points
    end

    min_y = points[1][2]
    min = 1
    
    # Bottom-most or else the left-most point when tie
    for i in 2:num_points
        y = points[i][2]
        if (y < min_y) || (y == min_y && points[i][1] < points[min][1])
            min_y = y
            min = i
        end
    end
    
    # Placing the bottom/left-most point at first position
    points[1], points[min] = points[min], points[1]

    global chull_p1 = points[1]

    points = Base.sort(points, lt = QCO._lt_filter)

    # If two or more points are collinear with chull_p1, remove all except the farthest from chull_p1
    arr_size = 2
    for i in 2:num_points
        while (i < num_points) && (QCO._orientation(chull_p1, points[i], points[i+1]) == 0)
            i += 1
        end
        points[arr_size] = points[i]
        arr_size += 1
    end

    if arr_size <= 3
        return []
    end

    chull = []
    append!(chull, [points[1]])
    append!(chull, [points[2]])
    append!(chull, [points[3]])
    
    for i in 4:(arr_size-1)
        while (size(chull)[1] > 1) && (QCO._orientation(chull[end - 1], chull[end], points[i]) != 2)
            pop!(chull)
        end
        append!(chull, [points[i]])
    end
    
    return chull
end

abstract AbstractMesh
abstract AbstractMesh1d <: AbstractMesh

typealias BCnames Vector{String}
"""
Stores an ordered 1D mesh
"""
type Mesh1d{T <: Real} <: AbstractMesh1d
    "Number of nodes"
    nb::Int
    "Number of elements"
    nel::Int
    "Coordinates of the nodes"
    x::Vector{T}
    "Index of neighbor elements"
    neigh::Array{Int,2}
    bcs::BCnames
end


function check_elem_order{T<:Real}(x::AbstractVector{T})
    n = length(x)

    for i = 2:n
        if x[i] <= x[i-1]
            return false
        end
    end
    return true
end

"""
Outer constructer for `Mesh1d` objects

It builds a mesh from the nodes
 
 * `x`: Vector with node positions
 * `bcs`: Vector with boundary condition tags

Alternatively, the mesh can be built from the end points 
and number of elements:

 * `a`: Left endpoint of the 1d domain
 * `b`: Right endpoint of the 1d domain
 * `nel`: Number of elements that the domain will be divided into
 * `bcs`: Vector with boundary condition tags
"""
function Mesh1d{T<:Real}(x::AbstractVector{T}, bcs=[])
    nb::Int = length(x)
    if nb < 2
        throw(DimensionMismatch("There should be at least 2 nodes in the 1d mesh"))
    end
    
    nel = nb-1

    if !check_elem_order(x)
        throw(ArgumentError("Nodes should be sequential from smaller to larger values"))
    end
                  
    y = [xx for xx in x]

    if length(bcs) > 2
        throw(DimensionMismatch("There should be at most one BC for each end of the domain"))
    end
    
    bcs2 = String[utf8(b) for b in bcs]

    neigh = zeros(Int, 2, nel)

    for e = 1:nel-1
        neigh[1,e] = e-1
        neigh[2,e] = e+1
    end
    neigh[1,nel] = nel-1
    neigh[2,nel] = 0

    return Mesh1d(nb, nel, y, neigh, bcs2)
    
end


function Mesh1d{T<:Real}(a::T, b::T, nel::Int, bcs=[])

    if a >= b
        throw(ArgumentError("Nodes should be sequential from smaller to larger values"))
    end
    
    nb = nel + 1

    x = linspace(a, b, nb)

    return Mesh1d(x, bcs)
    
end

"""
Stores information about an element.
"""
immutable mshElem1d{T<:Real}
    id::Int
    a::T
    b::T
    neigh::Tuple{Int,Int}
end

import Base.getindex
getindex(msh::Mesh1d, e::Int) = mshElem1d(e, msh.x[e], msh.x[e+1], (msh.neigh[1,e], msh.neigh[2,e]))


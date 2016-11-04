abstract AbstractCurve


type Quad{T<:Number}
    id::Int
    x::Vector{T}
    y::Vector{T}
    neigh::Vector{Tuple{Int,Int}}
end
nverts(::Type=Quad) = 4
nedges(::Type=Quad) = 4

nverts(el::Quad) = 4
nedges(el::Quad) = 4

xverts(el::Quad) = el.x
yverts(el::Quad) = el.y

neigh(el::Quad, e::Int) = el.neigh[e]


vertex_in_edge = [(1,2), (2,3), (3,4), (4,1)]
edge_in_vertex = [(4,1), (1,2), (2,3), (3,4)]

vinedge(e) = vertex_in_edge[e]
einvert(v) = edge_in_vertex[v]




type mshBndry{T<:Number}
    tag::String
    elem::Int
    face::Int
    params::Vector{T}
    funs::Vector{String}
end

type OtherInfo{T<:Number}
    name::String
    info::Vector{String}
end

type Mesh2d{T<:Number}
    elems::Vector{Quad{T}}
    bcs::Vector{mshBndry{T}}
    other::Vector{OtherInfo{T}}
end
    
    

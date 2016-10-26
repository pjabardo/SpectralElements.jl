using Jacobi

type Basis1d{T<:Number}
    "Number of nodes"
    Q::Int

    "Node coordinates in standard element"
    ξ::Vector{T}

    "Integration weights corresponding to each node"
    w::Vector{T}

    "Derivative matrix"
    D::Array{T,2}

    "Local numbering system"
    lnum::LocalNumSys1d
    
end

function Basis1d{T<:Number}(Q, ::Type{T}=Float64)
    ξ = zglj(Q, 0, 0, T)
    w = wglj(ξ)
    D = dglj(ξ)
    Basis1d{T}(Q, ξ, w, D, LocalNumSys1d(Q))
end

nquad(b::Basis1d) = b.Q
nmodes(b::Basis1d) = b.Q
qnodes(b::Basis1d) = b.ξ
qweights(b::Basis1d) = b.w
qdiffmat(b::Basis1d) = b.D
       

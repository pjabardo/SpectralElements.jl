using LinearAlgebra
using Jacobi

struct Basis1d{T<:Number}
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

function Basis1d(Q, ::Type{T}=Float64) where {T<:Number}
    ξ = zglj(Q, 0, 0, T)
    w = wglj(ξ)
    D = dglj(ξ)
    Basis1d{T}(Q, ξ, w, D, LocalNumSys1d(Q))
end

nquad(b::Basis1d) = b.Q
nmodes(b::Basis1d) = b.Q
qnodes(b::Basis1d) = b.ξ
qnodes(b::Basis1d, idx::Integer) = b.ξ[idx]
qweights(b::Basis1d) = b.w
qweights(b::Basis1d, idx::Integer) = b.w[idx]
qdiffmat(b::Basis1d) = b.D

lnum(b::Basis1d) = b.lnum


function ∂ξ!(b::Basis1d{T}, x::AbstractVector{T}, y::AbstractVector{T}) where {T<:Number}
    mul!(y, qdiffmat(b), x)
end
∂ξ(b::Basis1d{T}, x::AbstractVector{T}) where {T<:Number} = ∂ξ!(b, x, similar(x)) 

integrate(b::Basis1d{T}, x::AbstractVector{T}) where {T<:Number} =  dot(qweights(b), x)

        

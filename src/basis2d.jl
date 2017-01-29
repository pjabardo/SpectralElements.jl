

type Basis2d{T<:Number}
    "Number of nodes"
    Q::Int

    "Node coordinates in standard element"
    ξ::Vector{T}

    "Integration weights corresponding to each node"
    w::Vector{T}

    "Derivative matrix"
    D::Array{T,2}

    "Local numbering system"
    lnum::LocalNumSys2d
end


function Basis2d{T<:Number}(Q::Integer, ::Type{T}=FLoat64)

    ξ = zglj(Q, 0, 0, T)
    w = wglj(ξ)
    D = dglj(ξ)
    Basis2d{T}(Q, ξ, w, D, LocalNumSys2d(Q))
end
               


nquad(b::Basis2d) = b.Q
nmodes(b::Basis2d) = b.Q*b.Q
qnodes(b::Basis2d) = b.ξ
qweights(b::Basis2d) = b.w
qdiffmat(b::Basis2d) = b.D
       
function ∂ξ!{T<:Number}(b::Basis1d{T}, x::AbstractVector{T}, y::AbstractVector{T})
    A_mul_B!(y, qdiffmat(b), x)
end
∂ξ{T<:Number}(b::Basis1d{T}, x::AbstractVector{T}) = ∂ξ!(b, x, similar(x))

integrate{T<:Number}(b::Basis1d{T}, x::AbstractVector{T}) dot(qweights(b), x)

        

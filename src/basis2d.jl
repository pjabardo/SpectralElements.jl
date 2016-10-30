type Basis2d{T<:Number}
    "Number of modes in the first direction"
    Q::Tuple{Int,Int}
    "Node coordinates"
    ξ::Tuple{Vector{T},Vector{T}}
    "Integration weights"
    w::Tuple{Vector{T},Vector{T}}
    "Derivative matrix"
    D::Tuple{Matrix{T},Matrix{T}}
    "Local numbering system"
    #lnum::LocalNumSys2d
end

function Basis2d{T<:Number}(Q, ::Type{T}=Float64)
    ξ = zglj(Q, 0, 0, T)
    w = wglj(ξ)
    D = dglj(ξ)
    Basis2d{T}((Q, copy(Q)), (ξ, copy(ξ)),  (w, copy(w)), (D, copy(D)))#, LocalNumSys1d(Q))
end

function Basis2d{T<:Number}(Q₁, Q₂, ::Type{T}=Float64)
    ξ₁ = zglj(Q₁, 0, 0, T)
    w₁ = wglj(ξ₁)
    D₁ = dglj(ξ₁)

    ξ₂ = zglj(Q₂, 0, 0, T)
    w₂ = wglj(ξ₂)
    D₂ = dglj(ξ₂)

    Basis2d{T}( (Q₁, Q₂), (ξ₁, ξ₂), (w₁, w₂), (D₁, D₂) )
    
end

nquad(b::Basis2d) = b.Q
nquad(b::Basis2d, i) = b.Q[i]

nmodes(b::Basis2d) = b.Q
nmodes(b::Basis2d, i) = b.Q[i]

qnodes(b::Basis2d) = b.ξ
qnodes(b::Basis2d, i) = b.ξ[i]

qweights(b::Basis2d) = b.w
qweights(b::Basis2d, i) = b.w[i]

qdiffmat(b::Basis2d) = b.D
qdiffmat(b::Basis2d,i) = b.D[i]



    
    
    

   
    
    
    

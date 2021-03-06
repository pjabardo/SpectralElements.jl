abstract type AbstractElement end

struct Element1d{T<:Number} <: AbstractElement
    a::T
    b::T
    x::Vector{T}
    Jw::Vector{T}
    dξdx::Vector{T}
    bas::Basis1d{T}
end


function Element1d(a::T, b::T, bas::Basis1d{T}) where {T<:Number}
    w = qweights(bas)
    ξ = qnodes(bas)
    Q = nquad(bas)

    x = zeros(T, Q)
    dξdx = zeros(T,Q)
    Jw = zeros(T,Q)

    d = (b-a) / 2
        
    for i = 1:Q
        x[i] = (one(T)-ξ[i])*a/2 + (one(T) + ξ[i])*b/2
        dξdx[i] = one(T)/d
        Jw[i] = w[i] * d
    end
    
    Element1d(a, b, x, Jw, dξdx, bas)
end

basis(el::Element1d) = el.bas

nquad(el::Element1d) = nquad(el.bas)
nmodes(el::Element1d) = nquad(el.bas)
qnodes(el::Element1d) = qnodes(el.bas)
qweights(el::Element1d) = qweights(el.bas)
qdiffmat(el::Element1d) = qdiffmat(el.bas)




function ∂x!(el::Element1d, u, du)
    b = basis(el)
    D = qdiffmat(b)
    Q = nquad(b)
    dξ = el.dξdx
    for i = 1:Q
        for k = 1:Q
            du[k] += D[k,i] * u[i] * dξ[i]
        end
    end
    return du
end

∂x(el::Element1d, u) = ∂x!(el, u, similar(u))


function integrate(e::Element1d{T}, x::AbstractVector{T}) where {T<:Number}
    I = zero(T)
    
    N = nquad(e)
    Jw = e.Jw

    for i in 1:N
        I += Jw[i] * x[i]
    end

    return I
end

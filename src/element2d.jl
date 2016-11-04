

type Element2d{T<:Number} <: Element
    x::Array{T,2}
    y::Array{T,2}
    dξdx::Array{T,2}
    dξdy::Array{T,2}
    dηdx::Array{T,2}
    dηdy::Array{T,2}
    Jw::Array{T,2}
    Sjac::Vector{Vector{T}}
    vx::Vector{T}
    vy::Vector{T}
    is_curved::Bool
    bas::Basis1d{T}
    lnum::LocalNumSys2d
end


function Element2d{T<:Number}(bas::Basis1d{T}, lnum::LocalNumSys2d)

    Q = nquad(bas)

    x = zeros(T, Q, Q)
    y = zeros(T, Q, Q)
    dξdx = ones(T, Q, Q)
    dξdy = zeros(T, Q, Q)
    dηdx = zeros(T, Q, Q)
    dηdy = ones(T, Q, Q)
    Jw = zeros(T, Q, Q)
    Sjac = [zeros(T,Q) for e = 1:nedges(lnum)]
    vx = zeros(T, nverts(lnum))
    vy = zeros(T, nverts(lnum))

    Element2d{T}(x, y, dξdx, dξdy, dηdx, dηdy, Jw, Sjac, vx, vy, false, bas, lnum)

end

function Element2d{T<:Number}(bas::Basis1d{T}, lnum::LocalNumSys2d, vx, vy)
    el = Element2d{T}(bas, lnum)
    

end


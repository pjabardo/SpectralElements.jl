abstract Element

type Element1d{T<:Number} <: Element
    id::Int
    a::T
    b::T
    bas::Basis1d{T}
    J::Vector{T}
    wJ::Vector{T}
    dξdx::Vector{T}
end



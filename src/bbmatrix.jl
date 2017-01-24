"Abstract type that deals with global boundary-boundary matrices"
abstract BBSolver

type BBMatrix{T<:Number} <: BBSolver
    nb::Int
    nbslv::Int
    A::Array{T,2}
    fact::Base.LinAlg.LU{T,Array{T,2}}
    function BBMatrix(dof::DofMap)
        nb = nbmodes(dof)
        nbslv = nbslvmodes(dof)
        new(nb, nbslv, zeros(T, nbslv, nbslv))
    end
end

function assemble!{T<:Number}(Ag::BBMatrix{T}, Ae, m)

    nm = length(m)
    n1 = Ag.nbslv+1
    for i = 1:nm
        ig = m[i]
        if ig < n1
            for k = 1:nm
                kg = m[k]
                if kg < n1
                    Ag.A[kg,ig] += Ae[k,i]
                end
            end
        end
    end
end

function trf!(Ag::BBMatrix)
    Ag.fact = lufact!(Ag.A)
end
trs!(Ag::BBMatrix, x) = A_ldiv_B!(Ag.fact, x) 


"Global boundary-boundary Symmetric Positive-Definite matrix"
abstract BBSym <: BBSolver


type BBSymMatrix{T<:Number} <: BBSym
    nb::Int
    nbslv::Int
    A::Array{T,2}
    fact::Base.LinAlg.Cholesky{T,Array{T,2}}
    function BBSymMatrix(dof::DofMap)
        nb = nbmodes(dof)
        nbslv = nbslvmodes(dof)
        new(nb, nbslv, zeros(T, nbslv, nbslv))
    end

end

function assemble!{T<:Number}(Ag::BBSymMatrix{T}, Ae, m)

    nm = length(m)
    n1 = Ag.nbslv+1
    for i = 1:nm
        ig = m[i]
        if ig < n1
            for k = 1:nm
                kg = m[k]
                if kg < n1
                    Ag.A[kg,ig] += Ae[k,i]
                end
            end
        end
    end
end

function trf!(Ag::BBSymMatrix)
    Ag.fact = cholfact!(Ag.A)
end
trs!(Ag::BBSymMatrix, x) = A_ldiv_B!(Ag.fact, x) 

"Global boundary-boundary symmetric tridiagonal matrices"
type BBSymTri{T<:Number} <: BBSym
    "Number of boundary modes"
    nb::Int
    "Number of boundary modes that should be solved"
    nbslv::Int
    "Diagonal elements of the matrix"
    D::Array{T,1}
    "Sub-diagonal elements"
    E::Array{T,1}

    "Tridiagonal Matrix"
    tri::SymTridiagonal{T}
    "Factorization of symmetric tridiagonal matrix"
    fact:: LDLt{T,SymTridiagonal{T}}
    
    function BBSymTri(dof::DofMap)
        nb = nbmodes(dof)
        nbslv = nbslvmodes(dof)
        D = zeros(T,nbslv)
        E = zeros(T, max(nbslv-1,0))
        tri = SymTridiagonal(D, E)
        new(nb, nbslv, D, E, tri)
    end
end


"LU (or Choleksy) decomposition of boundary-boundary system"
function trf!(Ag::BBSymTri)
    Ag.fact = ldltfact!(Ag.tri)
end

"Solving linear system using LU (or Cholesky) decomposition"
trs!(Ag::BBSymTri, x) = A_ldiv_B!(Ag.fact, x)

"""
Assembles the global boundary-boundary matrix.
"""
function assemble!{T}(Ag::BBSymTri{T}, Ae, m)
    np1 = Ag.nbslv + 1
    D = Ag.D
    E = Ag.E

    ig = m[1]
    kg = m[2]

    if ig < np1
        D[ig] += Ae[1,1]
        if kg < np1
            D[kg] += Ae[2,2]
            E[min(ig,kg)] += Ae[1,2]
        end
    elseif kg < np1
        D[kg] += Ae[2,2]
    end

end



"Global boundary-boundary tridiagonal matrices"
type BBTri{T<:Number} <: BBSolver
    "Number of boundary modes"
    nb::Int
    "Number of boundary modes that should be solved"
    nbslv::Int
    "Diagonal elements of the matrix"
    D::Array{T,1}
    "Lower sub-diagonal elements"
    Dl::Array{T,1}
    "Upper sub-diagonal elements"
    Du::Array{T,1}
    "Tridiagonal object"
    tri::Tridiagonal{T}
    fact::Base.LinAlg.LU{T,Tridiagonal{T}}
    """
    Inner constructor. First the matrix must be assembled.
    """
    function BBTri(dof::DofMap)
        nb = nbmodes(dof)
        nbslv = nbslvmodes(dof)
        D = zeros(T,nbslv)
        Dl = zeros(T, max(nbslv-1, 0))
        Du = zeros(T, max(nbslv-1, 0))
        tri = Tridiagonal(Dl, D, Du)
        new(nb, nbslv, D, Dl, Du, tri)
    end
end


function assemble!{T}(Ag::BBTri{T}, Ae, m)
    np1 = Ag.nbslv + 1
    D = Ag.D
    Dl = Ag.Dl
    Du = Ag.Du


    ig = m[1]
    kg = m[2]

    if ig < np1
        D[ig] += Ae[1,1]
        if kg < np1
            D[kg] += Ae[2,2]
            if ig < kg
                Du[ig] += Ae[1,2]
                Dl[ig] += Ae[2,1]
            elseif kg < ig
                Du[kg] = Ae[2,1]
                Dl[kg] = Ae[1,2]
            end
        end
    elseif kg < np1
        D[kg] += Ae[2,2]
    end

end

function trf!(Ag::BBTri)
    Ag.fact = lufact!(Ag.tri)
    return
end
trs!(Ag::BBTri, x) = A_ldiv_B!(Ag.fact, x)




type BBTriP{T<:Number} <: BBSolver
    "Number of independent boundary modes"
    nb::Int
    "Number of boundary modes that should be solved"
    nbslv::Int
    "Lower diagonal"
    Dl::Vector{T}
    "Main Diagonal"
    D::Vector{T}
    "Upper diagonal"
    Du::Vector{T}
    "Auxiliary memory"
    x2::Vector{T}
    an1::T
    b1::T
    bn1::T
    cn::T
    cn1::T
    tri::Tridiagonal{T}
    fact::LinAlg.LU{T,Tridiagonal{T}}
    
    function BBTriP(dof::DofMap)
        nb = nbmodes(dof)
        nbslv = nbslvmodes(dof)

        nbslv = nb-1
        D = zeros(T,nbslv)
        Dl = zeros(T,nbslv-1)
        Du = zeros(T,nbslv-1)
        x2 = zeros(T,nbslv)
        tri = Tridiagonal(Dl, D, Du)
        new(nb, nbslv, Dl, D, Du, x2, zero(T), zero(T), zero(T), zero(T), zero(T), tri)
        
    end
            
end

function assemble!{T}(Ag::BBTriP{T}, Ae, m)

    m1 = m[1]
    m2 = m[2]
    
    nb = Ag.nb
    nbslv = Ag.nb-1

    ig = m1
    kg = m2

    D = Ag.D
    Dl = Ag.Dl
    Du = Ag.Du

    if m2 < nb && m2 != 1 #Not the last element 
        np1 = Ag.nbslv + 1
        
        
        if ig < np1
            D[ig] += Ae[1,1]
            if kg < np1
                D[kg] += Ae[2,2]
                if ig < kg
                    Du[ig] += Ae[1,2]
                    Dl[ig] += Ae[2,1]
                elseif kg < ig
                    Du[kg] = Ae[2,1]
                    Dl[kg] = Ae[1,2]
                end
            end
        elseif kg < np1
            D[kg] += Ae[2,2]
        end
    elseif m2 == nb # It is the next to last
        Ag.an1   += Ae[2,2]
        Ag.cn    += Ae[1,2]
        Ag.bn1   += Ae[2,1]
        D[nbslv] += Ae[1,1]
    elseif m2 == 1
        Ag.b1    += Ae[2,1]
        Ag.cn1   += Ae[1,2]
        D[1]     += Ae[2,2]
        Ag.an1   += Ae[1,1]
    end
    
end
                    

function trf!(Ag::BBTriP)
    Ag.fact = lufact!(Ag.tri)
    return
end
function trs!(Ag::BBTriP, x)
    n = Ag.nb
    n1 = n-1
    qn1 = view(x, 1:n1, :)
    nrhs = size(x,2)
    A_ldiv_B!(Ag.fact, qn1)
    x2 = Ag.x2
    for i = 1:nrhs
        fill!(x2, 0)
        x2[1] = -Ag.b1
        x2[n1] = -Ag.cn
        A_ldiv_B!(Ag.fact, x2)
        x[n,i] = (x[n,i] - Ag.cn1*qn1[1,i] - Ag.bn1*qn1[n1,i]) /
            (Ag.an1 + Ag.cn1*x2[1] + Ag.bn1*x2[n1])
        for k = 1:n1
            x[k,i] += x[n,i]*x2[k]
        end
    end
    x

end



type BBSymTriP{T<:Number} <: BBSym
    "Number of independent boundary modes"
    nb::Int
    "Number of boundary modes that should be solved"
    nbslv::Int
    "Main Diagonal"
    D::Vector{T}
    "Upper diagonal"
    Du::Vector{T}
    "Auxiliary memory"
    x2::Vector{T}
    an1::T
    cn::T
    cn1::T
    "Tridiagonal Matrix"
    tri::SymTridiagonal{T}
    "Factorization of symmetric tridiagonal matrix"
    fact:: LDLt{T,SymTridiagonal{T}}
    
    function BBSymTriP(dof::DofMap)
        nb = nbmodes(dof)
        nbslv = nbslvmodes(dof)

        nbslv = nb-1
        D = zeros(T,nb-1)
        Du = zeros(T,nb-2)
        x2 = zeros(T,nb-1)
        tri = SymTridiagonal(D, Du)
        new(nb, nbslv, D, Du, x2, zero(T), zero(T), zero(T), tri)
    end
end



function assemble!{T<:Number}(Ag::BBSymTriP{T}, Ae, m)

    m1 = m[1]
    m2 = m[2]
    
    nb = Ag.nb
    nbslv = Ag.nb-1

    ig = m1
    kg = m2

    D = Ag.D
    Du = Ag.Du
    
    if m2 < nb && m2 != 1 # Not the last element (if periodic)

        np1 = Ag.nbslv + 1
        
        
        if ig < np1
            D[ig] += Ae[1,1]
            if kg < np1
                D[kg] += Ae[2,2]
                if ig < kg
                    Du[ig] += Ae[1,2]
                elseif kg < ig
                    Du[kg] = Ae[2,1]
                end
            end
        elseif kg < np1
            D[kg] += Ae[2,2]
        end
    elseif m2 == nb # It is  the next to last
        Ag.an1   += Ae[2,2]
        Ag.cn    += Ae[1,2]
        D[nbslv] += Ae[1,1]
    elseif m2 == 1
        Ag.cn1   += Ae[1,2]
        D[1]     += Ae[2,2]
        Ag.an1   += Ae[1,1]
    end
    
end




function trf!(Ag::BBSymTriP)
    Ag.fact = ldltfact!(Ag.tri)
    return
end
function trs!(Ag::BBSymTriP, x)
    n = Ag.nb
    n1 = n-1
    qn1 = view(x, 1:n1, :)
    nrhs = size(x,2)
    A_ldiv_B!(Ag.fact, qn1)
    x2 = Ag.x2
    for i = 1:nrhs
        fill!(x2, 0)
        x2[1] = -Ag.cn1
        x2[n1] = -Ag.cn
        A_ldiv_B!(Ag.fact, x2)
        x[n,i] = (x[n,i] - Ag.cn1*qn1[1,i] - Ag.cn*qn1[n1,i]) /
            (Ag.an1 + Ag.cn1*x2[1] + Ag.cn*x2[n1])
        for k = 1:n1
            x[k,i] += x[n,i]*x2[k]
        end
    end
    x
end

@enum MatType FULLMAT TRIMAT TRIPMAT

"""
Generic matrix interface for 1d problems.

If the problem is very small, a full matrix will be used.
If it is larger and periodic, a BBTriP will be used. 
Else a BBTri is used.
"""
type BBMatrix1d{T<:Number} <: BBSolver
    nb::Int
    nbslv::Int
    mtype::MatType
    periodic::Bool
    mat::BBMatrix{T}
    tri::BBTri{T}
    trip::BBTriP{T}
    #BBMatrix1d(nb, nbslv, mtype, periodic) = new(nb, nbslv, mtype, periodic)
    function BBMatrix1d(dof::DofMap)
        nb = nbmodes(dof)
        nbslv = nbslvmodes(dof)
        periodic = false
        if dof.bmap[1,1] == dof.bmap[2,end]
            periodic = true
        end
        
        
        if nbslv <= 3
            mtype = FULLMAT
        elseif periodic == false
            mtype = TRIMAT
        else
            mtype = TRIPMAT
        end
        
        bb = new(nb, nbslv, mtype, periodic)

        if mtype == FULLMAT
            bb.mat = BBMatrix{T}(dof)
        elseif mtype == TRIMAT
            bb.tri = BBTri{T}(dof)
        else
            bb.trip = BBTriP{T}(dof)
        end
        
        return bb
    end
end


"""
Generic symmetric matrix interface for 1d problems.

If the problem is very small, a full matrix will be used.
If it is larger and periodic, a BBSymTriP will be used. 
Else a BBSymTri is used.
"""
type BBSymMatrix1d{T<:Number} <: BBSolver
    nb::Int
    nbslv::Int
    mtype::MatType
    periodic::Bool
    mat::BBSymMatrix{T}
    tri::BBSymTri{T}
    trip::BBSymTriP{T}
    function BBSymMatrix1d(dof::DofMap)
        nb = nbmodes(dof)
        nbslv = nbslvmodes(dof)
        periodic = false
        if dof.bmap[1,1] == dof.bmap[2,end]
            periodic = true
        end
        
        
        if nbslv <= 3
            mtype = FULLMAT
        elseif periodic == false
            mtype = TRIMAT
        else
            mtype = TRIPMAT
        end
        
        bb = new(nb, nbslv, mtype, periodic)

        if mtype == FULLMAT
            bb.mat = BBSymMatrix{T}(dof)
        elseif mtype == TRIMAT
            bb.tri = BBSymTri{T}(dof)
        else
            bb.trip = BBSymTriP{T}(dof)
        end
        
        return bb
    end
end



function assemble!{T<:Number}(Ag::Union{BBMatrix1d{T},BBSymMatrix1d{T}}, Ae, m)
    mt = Ag.mtype
    if mt == FULLMAT
        assemble!(Ag.mat, Ae, m)
    elseif mt == TRIMAT
        assemble!(Ag.tri, Ae, m)
    else
        assemble!(Ag.trip, Ae, m)
    end
    
end

function trf!(Ag::Union{BBMatrix1d,BBSymMatrix1d})
    mt = Ag.mtype
    if mt == FULLMAT
        trf!(Ag.mat)
    elseif mt == TRIMAT
        trf!(Ag.tri)
    else
        trf!(Ag.trip)
    end
end

function trs!(Ag::Union{BBMatrix1d,BBSymMatrix1d}, x) 
    mt = Ag.mtype
    if mt == FULLMAT
        trs!(Ag.mat,x)
    elseif mt == TRIMAT
        trs!(Ag.tri,x)
    else
        trs!(Ag.trip,x)
    end

end

import Base.LinAlg: BlasInt

type BBBandedMatrix{T<:Number} <: BBSolver
    nb::Int
    nbslv::Int
    bw::Int
    A::Matrix{T}
    ipiv::Vector{BlasInt}
    function BBBandedMatrix(dof::DofMap)

        bw = bandwidth(bw)
        nb = nbmodes(dof)
        nbslv = nbslvmodes(dof)
        nl = bw
        nu = bw
        A = zeros(T, 2*nl+nu+1, nbslv)
        ipiv = zeros(BlasInt,0)
        new(nb, nbslv, bw, A, ipiv)
    end        
end



function assemble!{T<:Number}(Ag::BBBandedMatrix{T}, Ae, m)

    nm = length(m)
    bw = Ag.bw
    n1 = Ag.nbslv+1
    A = Ag.A
    for i = 1:nm
        ig = m[i]
        if ig < n1
            for k = 1:nm
                kg = m[k]
                if kg < n1
                    A[bw+bw+1+kg-ig,ig] += Ae[k,i]
                end
            end
        end
    end
end


function trf!(Ag::BBBandedMatrix)

    A,Ag.ipiv = gbtrf!(Ag.bw, Ag.bw, Ag.nbslv, Ag.A)

end

trs!(Ag::BBBandedMatrix, x) = gbtrs!('N', Ag.bw, Ag.bw, Ag.nbslv, Ag.A, Ag.ipiv, x)


type BBSymBandedMatrix{T<:Number} <: BBSym
    nb::Int
    nbslv::Int
    bw::Int
    A::Matrix{T}
    function BBSymBandedMatrix(dof::DofMap)

        bw = bandwidth(bw)
        nb = nbmodes(dof)
        nbslv = nbslvmodes(dof)
        A = zeros(T, bw + 1, nbslv)
        new(nb, nbslv, bw, A)
    end        
end


function assemble!{T<:Number}(Ag::BBSymBandedMatrix{T}, Ae, m)

    nm = length(m)
    bw = Ag.bw
    n1 = Ag.nbslv+1
    A = Ag.A
    for i = 1:nm
        ig = m[i]
        if ig < n1
            for k = 1:nm
                kg = m[k]
                if kg < n1
                    if kg >= ig
                        A[1+kg-ig,ig] += Ae[k,i]
                    end
                end
            end
        end
    end
end

trf!(Ag::BBSymBandedMatrix) = pbtrf!('L', Ag.bw, Ag.A)
trs!(Ag::BBSymBandedMatrix, x) = pbtrs!('L', Ag.bw, Ag.A, x)

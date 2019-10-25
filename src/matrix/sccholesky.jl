

"Solver using Static condensation on symmetric positive definite problems."
mutable struct CholeskySC{T <: Number, Mat<:BBSolver, Dof <: DofMap} <: StaticCond
    dof::Dof
    Abb::Mat
    Aii::Vector{Matrix{T}}
    Aiifact::Vector{Base.LinAlg.Cholesky{T,Array{T,2}}}
    M::Vector{Matrix{T}}
    ub::Vector{T}
    Fie::Vector{Vector{T}}
    Fbe::Vector{Vector{T}}
    lft::Dict{Int,DirichiletLift}
    decomp::Bool
end

bbmatrix(solver::CholeskySC) = solver.Abb
dofmap(solver::CholeskySC) = solver.dof

function trf!(solver::StaticCond)
    if nbslvmodes(solver.dof) > 0
        trf!(solver.Abb)
    end
    solver.decomp = true
end



function add_local_matrix(solver::CholeskySC{T}, e::Integer, Abb::AbstractMatrix{T}, 
    Abi::AbstractMatrix{T},  Aib::AbstractMatrix{T}, Aii::AbstractMatrix{T}) where {T<:Number}

    dof = dofmap(solver)
    lmap = locmap(dof, e)
    nb = nbndry(lmap)
    ni = ninterior(lmap)
    ib = bndidx(lmap)

    Abb = zeros(T, nb, nb) 
    solver.Fbe[e] = zeros(T, nb)
    for i in 1:nb
        for k in 1:nb
            Abb[k,i] = Ae[ib[k], ib[i]]
        end
    end
        
    
    if ni > 0
        ii = intidx(lmap)
        Aii = zeros(T, ni, ni)
        Abi = zeros(T, nb, ni)
        M = zeros(T, ni, nb)
        solver.Fie[e] = zeros(T, ni)
        for i in 1:ni
            for k in 1:nb
                Abi[k,i] = Ae[ib[k], ii[i]]
            end
        end
        
        for i = 1:ni
            for k = 1:ni 
                Aii[k,i] = Ae[ ii[k], ii[i] ]
            end
        end
        fact = cholfact!(Aii)
        for k = 1:nb
            for i = 1:ni
                M[i,k] = Ae[ii[i],ib[k]]
            end
        end
        
        A_ldiv_B!(fact, M)

        solver.Aii[e] = Aii
        solver.Aiifact[e] = fact
        solver.M[e] = M
        
        gemm!('N', 'N', -one(T), Abi, M, one(T), Abb)

    end
    
    if hasdirbc(dof, e)
        solver.lft[e] = DirichiletLift(Abb, idirbc(dof, e))
    end
    
    assemble!(bbmatrix(solver), Abb, bmap(dof, e))
end


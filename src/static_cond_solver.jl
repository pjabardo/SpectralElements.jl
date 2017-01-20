"Abstract base class for linear system solvers"
abstract LinearSolver


"Abstract class for solvers using static condensation"
abstract StaticCond <: LinearSolver


"Solver using Static condensation on symmetric positive definite problems."
type CholeskySC{T <: Number, Mat<:BBSolver, Dof <: DofMap} <: StaticCond
    dof::Dof
    Abb::Mat
    Aii::Vector{Matrix{T}}
    Aiifact::Vector{Base.LinAlg.LU{T,Array{T,2}}}
    M::Vector{Matrix{T}}
    ub::Vector{T}
    Fi::Matrix{T}
    lft::Dict{Int,DirichiletLift}
    decomp::Bool
    Abbe::Matrix{T}
    Abie::Matrix{T}
    fbe::Vector{T}
end

bbmatrix(solver::CholeskySC) = solver.Abb
dofmap(solver::CholeskySC) = solver.dof

function trf!(solver::CholeskySC)
    if nbslvmodes(solver.dof) > 0
        trf!(solver.Abb)
    end
    solver.decomp = true
end



function CholeskySC{T<:Number, Mat<:BBSolver, Dof <: DofMap}(dof::Dof, ::Type{Mat}, 
                                                             ::Type{T}=Float64)
    nel = num_elems(dof)
    nbslv  = nbslvmodes(dof)
    nb = nbmodes(dof)
    Abb = Mat{T}(dof)
    Aii = Vector{Array{T,2}}(nel)
    #Aiifact = Vector{Base.LinAlg.Cholesky{T,Array{T,2}}}(nel)
    Aiifact = Vector{Base.LinAlg.LU{T,Array{T,2}}}(nel)
    M = Vector{Array{T,2}}(nel)

    lmap = locmap(dof)

    nbe = nbndry(lmap)
    nie = ninterior(lmap)
    Fi = zeros(T, nie, nel)
    for i = 1:nel
        Aii[i] = zeros(T, nie, nie)
        M[i] = zeros(T, nie, nbe)
    end
    ub = zeros(T, nbslv)
    lft = Dict{Int,DirichiletLift}()
    Abbe = zeros(T, nbe, nbe)
    Abie = zeros(T, nbe, nie)
    fbe = zeros(T, nbe)
    CholeskySC(dof, Abb, Aii, Aiifact, M, ub, Fi, lft, false, Abbe, Abie, fbe)
    
end

#using Base.LinAlg.BLAS.gemm!
#using Base.LinAlg.BLAS.gemv!
#using Base.LinAlg.LAPACK.potrf!
#using Base.LinAlg.LAPACK.potrs!

function add_local_matrix{T<:Number}(solver::CholeskySC{T}, e::Integer,
                                     Ae::AbstractMatrix{T})
    dof = dofmap(solver)
    lmap = locmap(dof)
    nb = nbndry(lmap)
    ni = ninterior(lmap)
    ib = bndidx(lmap)
    ii = intidx(lmap)

    Aii = solver.Aii[e]
    for i = 1:ni
        for k = 1:ni 
            Aii[k,i] = Ae[ ii[k], ii[i] ]
        end
    end
    #potrf!('L', Aii)
    #fact = cholfact!(Aii)
    fact = lufact!(Aii)
    solver.Aiifact[e] = fact
    M = solver.M[e]
    for k = 1:nb
        for i = 1:ni
            M[i,k] = Ae[ii[i],ib[k]]
        end
    end

    A_ldiv_B!(fact, M)
    Abb = solver.Abbe
    Abi = solver.Abie
    for i in 1:nb
        for k in 1:nb
            Abb[k,i] = Ae[ib[k], ib[i]]
        end
    end
    
    for i in 1:ni
        for k in 1:nb
            Abi[k,i] = Ae[ib[k], ii[i]]
        end
    end
            
    gemm!('N', 'N', -one(T), Abi, M, one(T), Abb)

    if hasdirbc(dof, e)
        solver.lft[e] = DirichiletLift(Abb, idirbc(dof, e))
    end
    
    assemble!(bbmatrix(solver), Abb, bmap(dof, e))
end


function solve!{T<:Number}(solver::CholeskySC{T}, Fe::AbstractMatrix{T})
    if !solver.decomp
        trf!(solver)
        solver.decomp = true
    end
    dof = dofmap(solver)
    lmap = locmap(dof)
    ib = bndidx(lmap)
    ii = intidx(lmap)
    nel = num_elems(dof)

    nbe = nbndry(lmap)
    nie = ninterior(lmap)
    Fb = solver.ub
    nb = nbmodes(dof)
    nbslv = nbslvmodes(dof)
    for i = 1:nbslv
        Fb[i] = zero(T)
    end
    Fbe = solver.fbe 
    for e = 1:nel

        Fie = view(solver.Fi, :, e)

        for i = 1:nbe
            Fbe[i] = Fe[ib[i],e]
        end
        if hasdirbc(dof, e)
            lift!(solver.lft[e], Fbe)
        end
        for i = 1:nie
            Fie[i] = Fe[ii[i],e]
        end
        gemv!('T', -one(T), solver.M[e], Fie, one(T), Fbe)
        
        m = bmap(dof, e)

        for i in 1:nbe
            ig = m[i]
            if ig <= nbslv
                Fb[ig] += Fbe[i]
            end
        end

    end


    # Solve linear system (boundary-boundary system
    Abb = bbmatrix(solver)
    if nbslv > 0
        trs!(Abb, Fb)
    end
    # Scatter the results and solve for each element:
    for e = 1:nel
        m = bmap(dof, e)
        for i = 1:nbe
            ig = m[i]
            if ig <= nbslv
                Fbe[i] = Fb[ig]
            else
                # Dirichilet BC
                Fbe[i] = Fe[ib[i],e]
            end
        end
        Fie = view(solver.Fi, :, e)
        #potrs!('L', solver.Aii[e], Fie)
        A_ldiv_B!(solver.Aiifact[e], Fie)
        gemv!('N', -one(T), solver.M[e], Fbe, one(T), Fie)

        for i = 1:nbe
            Fe[ib[i], e] = Fbe[i]
        end
        for i = 1:nie
            Fe[ii[i], e] = Fie[i]
        end
    end

    return Fe
    
end

               

               
type LU_SC{Mat<:BBSolver, T <: Number} <: StaticCond
    dof::DofMap
    Abb::Mat
    Aii::Array{T, 3}
    M::Array{T,3}
end







"Abstract base class for linear system solvers"
abstract LinearSolver


"Abstract class for solvers using static condensation"
abstract StaticCond <: LinearSolver


"Solver using Static condensation on symmetric positive definite problems."
type CholeskySC{T <: Number, Mat<:BBSolver, Dof <: DofMap} <: StaticCond
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

function trf!(solver::CholeskySC)
    if nbslvmodes(solver.dof) > 0
        trf!(solver.Abb)
    end
    solver.decomp = true
end


function add_local_matrix{T<:Number}(solver::CholeskySC{T}, e::Integer,
                                     Ae::AbstractMatrix{T})
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


function solve!{T<:Number}(solver::CholeskySC{T}, Fe::AbstractMatrix{T})
    if !solver.decomp
        trf!(solver)
        solver.decomp = true
    end
    dof = dofmap(solver)
    nel = num_elems(dof)


    Fb = solver.ub
    nb = nbmodes(dof)
    nbslv = nbslvmodes(dof)


    for i = 1:nbslv
        Fb[i] = zero(T)
    end
    for e = 1:nel

        lmap = locmap(dof, e)
        nbe = nbndry(lmap)
        nie = ninterior(lmap)

        ib = bndidx(lmap)
        ii = intidx(lmap)

        Fie = solver.Fie[e]
        Fbe = solver.Fbe[e]

        for i = 1:nbe
            Fbe[i] = Fe[ib[i],e]
        end
        if hasdirbc(dof, e)
            lift!(solver.lft[e], Fbe)
        end

        if nie > 0
            for i = 1:nie
                Fie[i] = Fe[ii[i],e]
            end
            gemv!('T', -one(T), solver.M[e], Fie, one(T), Fbe)
        end
        
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
        Fbe = solver.Fbe[e]
        for i = 1:nbe
            ig = m[i]
            if ig <= nbslv
                Fbe[i] = Fb[ig]
            else
                # Dirichilet BC
                Fbe[i] = Fe[ib[i],e]
            end
        end

        if nie > 0
            Fie = solver.Fie[e]
            A_ldiv_B!(solver.Aiifact[e], Fie)
            gemv!('N', -one(T), solver.M[e], Fbe, one(T), Fie)
        end
        
        for i = 1:nbe
            Fe[ib[i], e] = Fbe[i]
        end

        if nie > 0
            for i = 1:nie
                Fe[ii[i], e] = Fie[i]
            end
        end
    end

    return Fe
    
end

               


               
type LUSC{T<:Number, Mat<:BBSolver, Dof<:DofMap} <: StaticCond
    dof::Dof
    Abb::Mat
    Aii::Vector{Matrix{T}}
    Aiifact::Vector{Base.LinAlg.LU{T,Array{T,2}}}
    M1::Vector{Matrix{T}}
    M2::Vector{Matrix{T}}
    ub::Vector{T}
    Fie::Vector{Vector{T}}
    Fbe::Vector{Vector{T}}
    lft::Dict{Int,DirichiletLift}
    decomp::Bool
end



function LUSC{T<:Number, Mat<:BBSolver, Dof <: DofMap}(dof::Dof, ::Type{Mat}, 
                                                       ::Type{T}=Float64)
    nel = num_elems(dof)
    nbslv  = nbslvmodes(dof)
    nb = nbmodes(dof)
    Abb = Mat{T}(dof)
    
    Aii = Vector{Array{T,2}}(nel)
    Fie = Vector{Vector{T}}(nel)
    Fbe = Vector{Vector{T}}(nel)

    Aiifact = Vector{Base.LinAlg.Cholesky{T,Array{T,2}}}(nel)
    M = Vector{Array{T,2}}(nel)

    ub = zeros(T, nbslv)

    lft = Dict{Int,DirichiletLift}()

    CholeskySC(dof, Abb, Aii, Aiifact, M, ub, Fie, Fbe, lft, false)
    
end

#using Base.LinAlg.BLAS.gemm!
#using Base.LinAlg.BLAS.gemv!
#using Base.LinAlg.LAPACK.potrf!
#using Base.LinAlg.LAPACK.potrs!

function add_local_matrix{T<:Number}(solver::LUSC{T}, e::Integer,
                                     Ae::AbstractMatrix{T})
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


function solve!{T<:Number}(solver::LUSC{T}, Fe::AbstractMatrix{T})
    if !solver.decomp
        trf!(solver)
        solver.decomp = true
    end
    dof = dofmap(solver)
    nel = num_elems(dof)


    Fb = solver.ub
    nb = nbmodes(dof)
    nbslv = nbslvmodes(dof)


    for i = 1:nbslv
        Fb[i] = zero(T)
    end
    for e = 1:nel

        lmap = locmap(dof, e)
        nbe = nbndry(lmap)
        nie = ninterior(lmap)

        ib = bndidx(lmap)
        ii = intidx(lmap)

        Fie = solver.Fie[e]
        Fbe = solver.Fbe[e]

        for i = 1:nbe
            Fbe[i] = Fe[ib[i],e]
        end
        if hasdirbc(dof, e)
            lift!(solver.lft[e], Fbe)
        end

        if nie > 0
            for i = 1:nie
                Fie[i] = Fe[ii[i],e]
            end
            gemv!('T', -one(T), solver.M[e], Fie, one(T), Fbe)
        end
        
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
        Fbe = solver.Fbe[e]
        for i = 1:nbe
            ig = m[i]
            if ig <= nbslv
                Fbe[i] = Fb[ig]
            else
                # Dirichilet BC
                Fbe[i] = Fe[ib[i],e]
            end
        end

        if nie > 0
            Fie = solver.Fie[e]
            A_ldiv_B!(solver.Aiifact[e], Fie)
            gemv!('N', -one(T), solver.M[e], Fbe, one(T), Fie)
        end
        
        for i = 1:nbe
            Fe[ib[i], e] = Fbe[i]
        end

        if nie > 0
            for i = 1:nie
                Fe[ii[i], e] = Fie[i]
            end
        end
    end

    return Fe
    
end

               









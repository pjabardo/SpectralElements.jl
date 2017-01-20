# DofMap
abstract DofMap


"""
Degree of freedom map for 1d discretizations

This type assumes that all elements have the same internal
discretization and same basis function. It also assumes that the 
basis function can be split into:
 
 * interior modes that are zero on the boundary
 * boundary modes that are nonzero at a sinlge boundary node.
"""
type DofMap1d <: DofMap
    "Number of boundary modes"
    nb::Int

    "Number of boundary modes that are not on Dirichilet BCs"
    nbslv::Int

    "Number of Dirichilet boundary modes"
    nbdir::Int

    "Number of elements"
    nel::Int

    "Local numbering system map for each element (the same for all elements)"
    lmap::LocalNumSys1d

    "Boundary map. Maps local boundary modes into global modes"
    bmap::Array{Int,2}

    "Maps local modes into global modes"
    map::Array{Int,2}

    "Identifies Dirichilet modes"
    idir::Dict{Int,Vector{Int}}

    "Inner constructor. Should not usually be called"
    function DofMap1d(nel, nb, nbdir, lmap::LocalNumSys1d,
                      bmap::Array{Int,2}, idir::Dict{Int,Vector{Int}})
        nbslv = nb - nbdir
        nloc = nmodes(lmap)
        mp = zeros(Int, nloc, nel)
        for e = 1:nel
            mp[1,e] = bmap[1,e]
            mp[2,e] = bmap[2,e]
            for i = 1:ninterior(lmap)
                mp[2+i,e] = nb + (e-1)*ninterior(lmap) + i
            end
        end
        new(nb, nbslv, nbdir, nel, lmap, bmap, mp, idir)
    end
end

"Number of boundary modes in the discretization"
nbmodes(dof::DofMap) = dof.nb
"Number of boundary modes in the discretization"
nbndry(dof::DofMap) = dof.nb

"Number of boundary modes that are not Dirichilet BCs"
nbslvmodes(dof::DofMap) = dof.nbslv
"Number of elements"
num_elems(dof::DofMap) = dof.nel

"Total number of interior modes"
ninodes(dof::DofMap1d) = dof.nel * ninterior(dof.lmap)
"Total number of interior modes"
ninterior(dof::DofMap1d) = dof.nel * ninterior(dof.lmap)

"Total number of modes not in Dirichilet BCs (interior and boundary modes)"
nslvmodes(dof::DofMap1d) = nbslvmodes(dof) + ninodes(dof)

"Returns the local numbering system"
locmap(dof::DofMap1d, e=1) = dof.lmap

"Is the element located on a Dirichilet boundary condition"
hasdirbc(dof::DofMap1d, e) = haskey(dof.idir, e)

"Returns the Dirichilet modes for a given element"
idirbc(dof::DofMap1d, e) = dof.idir[e]

export DofMap1d
export nbmodes, nbslvmodes, num_elems, ninodes, nslvmodes, locmap, hasdirbc, idirbc

"""
Local to global degree of freedom mapping

This function creates a `DofMap1d` object given a local mapping scheme, a the
number of nodes and the index of nodes on Dirichilet boundary conditions. It
also allows for the case where the boundary condition is periodic, such that the first node
equals to the last node.
"""
function DofMap1d(lmap, nnodes, idir, iper=false)
    nbe = nbndry(lmap)
    nie = ninterior(lmap)
    
    nel = nnodes - 1
    bmap = zeros(Int, nbe, nel)
    
    for e = 1:nel
        bmap[1,e] = e
        bmap[2,e] = e+1
    end
    nb = nnodes
    nd = length(idir)
    if iper
        bmap[2,nel] = 1
        nb = nnodes - 1
        nd = 0
    elseif nd > 0 # There is a Dirichilet BC.
        if nd == 1
            if idir[1] == 1
                ii = [nnodes;  1:(nnodes-1);]
            else
                ii = [1:nnodes;]
            end
        else
            ii = [nnodes; 1:(nnodes-1);]
        end
        
        for e = 1:nel
            bmap[1,e] = ii[e]
            bmap[2,e] = ii[e+1]
        end
        
    end
    idir = Dict{Int,Vector{Int}}()
    nbslv = nb - nd

    ib = bndidx(lmap)
    for e = 1:nel
        if bmap[1,e] > nbslv && bmap[2,e] > nbslv
            idir[e] = [1,2] 
        elseif bmap[1,e] > nbslv
            idir[e] =  [1] 
        elseif bmap[2,e] > nbslv
            idir[e] = [2] 
        end
    end
    return DofMap1d(nel, nb, nd, lmap, bmap, idir)
end

num_be(dof::DofMap, e) = nbndry(dof.lmap)
num_ie(dof::DofMap, e) = ninterior(dof.lmap)

bmap(dof, e) = view(dof.bmap, :, e)

function global2local(dof::DofMap, xg::Array{Float64,1})

    xloc = zeros(dof.nloc, dof.nel)
    for e = 1:dof.nel
        xloc[:,e] = xg[dof.map[:,e]]
    end
    return xloc
end


function bandwidth(dof::DofMap)

    bm = bmap(dof)
    nel = size(bm,2)
    nbe = size(bm,1)
    nb = nbmodes(dof)
    nbslv = nbslvmodes(dof)
    bw = 0
    min_idx = 0
    max_idx = 0
    for e=1:nel
        for i = 1:nbmodes
            b = bm[i,e]
            if b <= nbslv
                min_idx = b
                max_idx = b
                break
            end
        end
        for i=1:nbmodes
            b = bm[i,e]
            if b < nbslv
                min_idx = min(b, min_idx)
                max_idx = max(b, max_idx)
            end
            bw = max(bw, (max_idx-min_idx))
        end
    end
    return bw
end


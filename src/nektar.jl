


function section(flines, header)

    nl = length(flines)
    for i in 1:nl
        if length(search(flines[i], header)) > 0
            return i
        end
    end

    return -1
end
    

function element(rea, idx)
    if length(search(uppercase(rea[1]), "QUAD")) > 0
        el = Quad()
    else
        throw(ArgumentException("Only Quad type elements are possible"))
    end
    xverts = [parse(Float64, x) for x = split(strip(rea[2]))][1:nverts(el)]
    yverts = [parse(Float64, x) for x = split(strip(rea[2]))][1:nverts(el)]
    neigh = [(0,0) for i in 1:nedges(Quad)]

    Quad{Float64}(idx, xverts, yverts, neigh)
    
end
    
    
function mesh(rea)

    isec = section(rea, "MESH DATA") + 1
    nelem = parse(Int, split(strip(rea[isec]), " ")[1])
    elems = Array(Quad{Float64}, nelem)
    isec += 1
    for e = 1:nelem
        elems[e] = element(rea[isec:(isec+2)], e)
        isec += 3
    end

    return elems
    
end
        


function bndry_sec(rea, isec, nterms=2)

    s = split(strip(rea[isec]))
    tag = strip(s[1])
    nums = [parse(Float64, strip(ss)) for ss in s[2:end]]
    el = Int(round(nums[1]))
    ed = Int(round(nums[2]))
    p = nums[3:end]
    if islower(tag)
        funs = rea[(isec+1):(isec+nterms)]
        isec += nterms + 1
    else
        funs = String[]
        isec += 1
    end
    return mshBndry(tag, el, ed, p, funs), isec
end

function bndry_conds(rea, elems)
    isec = section(rea, "BOUNDARY CONDITIONS")+2
    bndrylst = mshBndry{Float64}[]
    for el in elems
        ne = Mesh2d.nedges(el)
        for edge = 1:ne
            b, isec = bndry_sec(rea, isec)
            if b.tag == 'E' || b.tag == 'P'
                el.neigh[1] = ( Int(round(b.params[1])), Int(round(b.params[2])) )
            else
                push!(bndrylst, b)
            end
        end
        if ne == 3
            isec += 1
        end
    end
    nb = length(bndrylst)
    for i = 1:nb
        el = bndrylst[i].elem
        ed = bndrylst[i].edge
        elems[el].neigh[ed] =(-1, -i)
    end
    return bndrylst
end
    

type NekCurve 
    cname::String
    line::String
end

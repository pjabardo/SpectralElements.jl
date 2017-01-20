
"""
Implements Dirichilet lift

In FEM, Neuman boundary conditions are naturally satisfied
(hey, they are also known as natural BCs!) but Dirichilet BCs
are another problem. The simplest solution is to use basis functions
that satisfy the boundary conditions. Now, by design, your solution
satisfies the BCs. This usually is a difficult thing to do. Except
when the solution is zero at Dirichilet BCs. This homogeneous case
is easy with basis that split into boundary modes and interior modes:

Interior modes satisy the homogeneous requirements! Just drop the relevant
boundary modes!

Now we still have the issue of non-homogeneous BCs. That is easy as well,
just split the solution into two parts: homogeneous and dirichilet:


u(x) = u^h(x) + u^d(x)


The weak form will take the following form:


a(u, v) = a(u^h + u^d, v) = a(u^h, v) + a(u^d, v) = f(x)


Since u^d is known (Dirichiulet BCs),


a(u^h, v) = f(x) - a(u^d, v)


The `DirichiletLift` type and associated methods implement the correction of
the right hand side of the equation at a local element.

This implementation is very simple and probably not the most efficient. Just
compute the matrix as if it were a common internal element. Store this matrix
and the indicies that correspond to the Dirichilet modes.

When you want to lift the solution, just set the correct Dirichilet modes
and the `lift!` function will set the other modes to zero and subtract
the contribution.

 The DirichiletLift constructor accepts these arguments

 * `A0` local elemental matrix
 * `idir` local index of dirichilet modes

"""
type DirichiletLift{T <: Number}
    "Order of local matrix"
    M::Int

    "Number of dof of non-Dirichilet modes"
    nh::Int

    "Number of dof of Dirichilet modes"
    nd::Int
    
    "Operator matrix (copy)"
    A::Array{T,2}

    "Index of non-Dirichilet modes"
    ih::Vector{Int}

    "Index of Dirichilet modes"
    id::Vector{Int}
    

end



"""
    Create a new DirichiletLift object.
    
    * `A0` Full elemental matrix
    * `idir` Index of Dirichilet modes
    """
function DirichiletLift{T<:Number}(A0::AbstractArray{T,2}, idir)
    M = size(A0,1)
    nd = length(idir)
    nh = M - nd
    A = zeros(T, nh, nd)
    ih = zeros(Int, nh)
    id = zeros(Int, nd)
    isx = falses(M)

    for i in 1:nd
        isx[idir[i]] = true
        id[i] = idir[i]
    end
    count = 1
    for i = 1:M
        if !isx[i]
            ih[count] = i
            count = count + 1
        end
    end
    for j = 1:nd
        for k = 1:nh
            A[k,j] = A0[ih[k], id[j]]
        end
    end
    
    DirichiletLift(M, nh, nd, A, ih, id)
end


"""
Actually lifts the solution

 * `lft` A `DirichiletLift` object that was previously created.
 * `f` A vector of the RHS. Dirichilet BCs are stored in it

"""
function lift!{T<: Number}(lft::DirichiletLift, f::AbstractVector{T})
    nd = lft.nd
    nh = lft.nh
    id = lft.id
    ih = lft.ih
    A = lft.A
    for i = 1:nd
        for k = 1:nh
            f[ih[k]] -= A[k,i] * f[id[i]]
        end
    end
end




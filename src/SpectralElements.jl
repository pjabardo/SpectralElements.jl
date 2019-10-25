module SpectralElements

# package code goes here

module Mesh
include("mesh1d.jl")
include("mesh2d.jl")
include("nektar.jl")

end

include("locnum.jl")
#include("basis1d.jl")
#include("element1d.jl")
#include("element2d.jl")

end # module

SEM = SpectralElements

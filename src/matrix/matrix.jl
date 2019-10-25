

#module Matrix

"Abstract base class for linear system solvers"
abstract type LinearSolver end


"Abstract class for solvers using static condensation"
abstract type StaticCond <: LinearSolver end


include("bbmatrix.jl")

#include("sccholesky.jl")

#end

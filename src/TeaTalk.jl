module TeaTalk

using Reexport

@reexport using LinearAlgebra
@reexport using CompScienceMeshes
@reexport using BEAST

@reexport using Plots

Base.ones(::Type{T}, m::Array) where {T} = Base.ones(T, size(m))

center = CompScienceMeshes.center

include("primitives.jl")

export center
export meshdisk

end # module

module TeaTalk

using Reexport

@reexport using LinearAlgebra
@reexport using CompScienceMeshes
@reexport using BEAST
@reexport using Plots

Base.ones(::Type{T}, m::Matrix) where {T} = Base.ones(T, size(m))

end # module

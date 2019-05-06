module TeaTalk

using Reexport

@reexport using LinearAlgebra
@reexport using CompScienceMeshes
@reexport using BEAST

ENV["PLOTS_DEFAULT_BACKEND"] = "plotlyjs"
@reexport using Plots

Base.ones(::Type{T}, m::Array) where {T} = Base.ones(T, size(m))

end # module

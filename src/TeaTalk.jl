module TeaTalk

using Reexport

@reexport using LinearAlgebra
@reexport using CompScienceMeshes
@reexport using BEAST
@reexport using Makie

Base.ones(::Type{T}, m::Array) where {T} = Base.ones(T, size(m))

center(mesh,cell) = cartesian(CompScienceMeshes.center(chart(mesh,cell)))

include("primitives.jl")

function expand_VF(m)
    @assert dimension(m) == 2
    @assert universedimension(m) == 3
    v = vertexarray(m)
    V = zeros(eltype(v), (3*numcells(m), 3))
    for (i,c) in enumerate(cells(m))
        @assert 1 <= c[1] <= numvertices(m)
        @assert 1 <= c[2] <= numvertices(m)
        @assert 1 <= c[3] <= numvertices(m)
        V[3(i-1)+1,:] = v[c[1],:]
        V[3(i-1)+2,:] = v[c[2],:]
        V[3(i-1)+3,:] = v[c[3],:]
    end
    F = [3(i-1)+j for i in 1:numcells(m), j in 1:3]
    return V,F
end

function Makie.convert_arguments(P::Type{<:Makie.Mesh}, m::CompScienceMeshes.Mesh)
    V,F = expand_VF(m)
    return Makie.convert_arguments(P,V,F)
end

fc(c) = (C=similar(c,3length(c)); for (i,x) ∈ pairs(c) C[3(i-1).+(1:3)].=x end; C)

export fc
export center
export meshdisk

end # module

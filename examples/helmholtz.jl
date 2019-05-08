using TeaTalk

Id = TeaTalk.Identity()

h = 0.05
m = meshdisk(radius=1.0, h=h)
# b = boundary(m)
#
# X = lagrangec0d1(m)
# Y = curl(X)
# A = assemble(Id,Y,Y)
# F = eigen(A)
# p = sortperm(F.values)
#
# u = F.vectors[:,p[3]]


# using PlotlyJS
# fcr, geo = facecurrents(u,X)
#
# PlotlyJS.plot(patch(m,getindex.(fcr,1)))

import Makie


function expand_VF(m)
    @show size(color)
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
    # FC = zeros(3*numcells(m))
    # for i in 1:numcells(m)
    #     FC[3(i-1)+1] = color[i]
    #     FC[3(i-1)+2] = color[i]
    #     FC[3(i-1)+3] = color[i]
    # end
    # @show size(FC)
    # @assert length(FC) == size(V,1)
    return V,F
end

function Makie.convert_arguments(P::Type{<:Makie.Mesh}, m::CompScienceMeshes.Mesh)
    @info "Converting args..."
    # @show size(color)
    # @assert dimension(m) == 2
    # @assert universedimension(m) == 3
    # v = vertexarray(m)
    # V = zeros(eltype(v), (3*numcells(m), 3))
    # for (i,c) in enumerate(cells(m))
    #     @assert 1 <= c[1] <= numvertices(m)
    #     @assert 1 <= c[2] <= numvertices(m)
    #     @assert 1 <= c[3] <= numvertices(m)
    #     V[3(i-1)+1,:] = v[c[1],:]
    #     V[3(i-1)+2,:] = v[c[2],:]
    #     V[3(i-1)+3,:] = v[c[3],:]
    # end
    # F = [3(i-1)+j for i in 1:numcells(m), j in 1:3]
    # FC = zeros(3*numcells(m))
    # for i in 1:numcells(m)
    #     FC[3(i-1)+1] = color[i]
    #     FC[3(i-1)+2] = color[i]
    #     FC[3(i-1)+3] = color[i]
    # end
    # @show size(FC)
    # @assert length(FC) == size(V,1)
    V,F = expand_VF(m)
    return Makie.convert_arguments(P,V,F)
end

coordinates = [
    0.0 0.0;
    0.5 0.0;
    1.0 0.0;
    0.0 0.5;
    0.5 0.5;
    1.0 0.5;
    0.0 1.0;
    0.5 1.0;
    1.0 1.0;]
connectivity = [
    1 2 5;
    1 4 5;
    2 3 6;
    2 5 6;
    4 5 8;
    4 7 8;
    5 6 9;
    5 8 9;]
color = [0.0, 0.0, 0.0, 0.0, -0.375, 0.0, 0.0, 0.0, 0.0]
scene = Makie.mesh(coordinates, connectivity, color = color, shading = false)

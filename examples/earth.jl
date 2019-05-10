using TeaTalk
d = pathof(TeaTalk)
Id = TeaTalk.Identity()

m = readmesh(joinpath(dirname(d),"../examples/meshes","sphere.in"))
m = meshsphere(1.0, 0.05)
@show numcells(m)
@show numvertices(m)

v = skeleton(m,0)
function ocean(cell)
    p = center(v,cell)
    p0 = point(1,0,0)
    r0 = 0.5
    return norm(p-p0) > r0
    # tol = eps(eltype(p))
    # abs(p[1] - 0.5) < tol && return false
    # abs(p[1] + 0.5) < tol && return false
    # abs(p[2] - 0.5) < tol && return false
    # abs(p[2] + 0.5) < tol && return false
    # p[1] > tol || p[2] > tol
end
v1 = submesh(ocean, v)



X = lagrangec0d1(m, getindex.(cells(v1),1), Val{3})
@show numfunctions(X)

function f(x)
    x0 = point(-1.0,1.0,1.0)
    x0 /= norm(x0)
    r = norm(x-x0)
    r0 = 0.1
    return exp(-(r/r0)^2)/sqrt(r0)
end

g = TeaTalk.BEAST.ScalarTrace(f)
κ = 4π
@show 2π/κ

@hilbertspace u
@hilbertspace v
Eq = @varform Id[curl(v),curl(u)] - κ^2*Id[v,u] ==g[v] #- im*κ*Id[trc(v),trc(u)] == g[v]
eq = @discretise Eq u∈X v∈X

uh = solve(eq)

fcr, geo = facecurrents(uh,X)
# c = clamp.(real.(getindex.(fcr,1)),-0.02,+0.02)
c = real.(getindex.(fcr,1))
# c = log10.(abs.(real.(getindex.(fcr,1))))
Makie.mesh(m, color=fc(c), shading=false)

using TeaTalk

Id = TeaTalk.Identity()

h = 0.05
m = meshdisk(radius=1.0, h=h)

# X = lagrangec0d1(m, getindex.(cells(v),1), Val{3})
X = lagrangec0d1(m)
Y = curl(X)
A = assemble(Id,Y,Y)
F = eigen(A)
p = sortperm(F.values)

u = F.vectors[:,p[4]]
fcr, geo = facecurrents(u,X)
Makie.mesh(m, color=fc(getindex.(fcr,1)), colormap=:lighttest, shading=:false)


h = 0.025
m = meshdisk(radius=1.0, h=h)

v = skeleton(m,0)
b = boundary(m)
X = lagrangec0d1(m, getindex.(cells(v),1), Val{3})
@show numfunctions(X)
x = strace(X,b)

function f(x)
    r = norm(x-point(0.3,0.3,0.0))
    r0 = 0.1
    return exp(-(r/r0)^2)/sqrt(r0)
end

g = TeaTalk.BEAST.ScalarTrace(f)
κ = 8π
trc = X->strace(X,b)

@hilbertspace u
@hilbertspace v
Eq = @varform Id[curl(v),curl(u)] - κ^2*Id[v,u] - im*κ*Id[trc(v),trc(u)] == g[v]
eq = @discretise Eq u∈X v∈X

uh = solve(eq)

fcr, geo = facecurrents(uh,X)
c = clamp.(real.(getindex.(fcr,1)),-0.002,+0.002)
# c = log10.(abs.(real.(getindex.(fcr,1))))
Makie.mesh(m, color=fc(c), shading=false)

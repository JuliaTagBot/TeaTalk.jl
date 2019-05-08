using TeaTalk

Id = TeaTalk.Identity()

h = 0.05
m = meshdisk(radius=1.0, h=h)
v = skeleton(m,0)
b = boundary(m)

# X = lagrangec0d1(m, getindex.(cells(v),1), Val{3})
X = lagrangec0d1(m)
Y = curl(X)
A = assemble(Id,Y,Y)
F = eigen(A)
p = sortperm(F.values)

u = F.vectors[:,p[4]]
fcr, geo = facecurrents(u,X)
Makie.mesh(m, color=fc(getindex.(fcr,1)), colormap=:lighttest, shading=:false)

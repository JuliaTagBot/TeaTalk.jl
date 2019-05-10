using TeaTalk

h = 0.025
m = meshrectangle(1.0, 1.0, h, 3)
nx = ny = round(Int,1/h-1)

X = lagrangec0d1(m)
Id = BEAST.Identity()

Y = curl(X)

A = assemble(Id,Y,Y)
F = eigen(A)
p = sortperm(F.values)
u = F.vectors[:,p[6]]

heatmap(reshape(u,(nx,ny)))
fcr, geo = facecurrents(u,X)
Makie.mesh(m, color=fc(getindex.(fcr,1)), shading=false)

# Try to recreate the matlab logo
using TeaTalk
Id = BEAST.Identity()
h = 0.025
m = meshrectangle(1.0, 1.0, h, 3)
m = CompScienceMeshes.translate(m, [-0.5, -0.5, 0.0])
nx = ny = round(Int,1/h-1)

verts = skeleton(m,0)
function notll(cell)
    p = center(verts,cell)
    tol = eps(eltype(p))
    abs(p[1] - 0.5) < tol && return false
    abs(p[1] + 0.5) < tol && return false
    abs(p[2] - 0.5) < tol && return false
    abs(p[2] + 0.5) < tol && return false
    p[1] > tol || p[2] > tol
end

verts1 = submesh(notll, verts)
X = lagrangec0d1(m, getindex.(cells(verts1),1), Val{3})
Y = curl(X)

A = assemble(Id,Y,Y)
F = eigen(A)
p = sortperm(F.values)
u = F.vectors[:,p[1]]

fcr, geo = facecurrents(u,X)
Makie.mesh(m, color=fc(getindex.(fcr,1)), shading=false)

H = zeros(nx+2,ny+2)
for (i,c) in pairs(cells(verts1))
    x = center(verts1,c)
    p = round(Int,(x[2]+0.5)/h)+1
    q = round(Int,(x[1]+0.5)/h)+1
    H[p,q] = u[i]
end
surface(H/maximum(H)*size(H,1)/2, shading=false)

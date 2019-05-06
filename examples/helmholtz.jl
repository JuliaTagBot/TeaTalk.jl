using TeaTalk

h = 0.025
m = meshrectangle(1.0, 1.0, h, 3)
nx = ny = round(Int,1/h-1)

X = lagrangec0d1(m)
Id = BEAST.Identity()

Y = curl(X)

A = assemble(Id,Y,Y)
F = eigen(A)

heatmap(reshape(real.(F.vectors[:,2]),(nx,ny)))

fcr, geo = facecurrents(ones(numfunctions(Y)),Y)
patch(m,norm.(fcr))

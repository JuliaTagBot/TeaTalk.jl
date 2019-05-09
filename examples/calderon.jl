using TeaTalk
Γ = readmesh(joinpath(dirname(d),"../examples/meshes","sphere2.in"))

X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

κ = ω = 1.0; γ = κ*im
N = BEAST.NCross()
T = Maxwell3D.singlelayer(wavenumber=κ)

Txx = assemble(T,X,X)
F_Txx = eigen(Txx)
@show maximum(abs.(F_Txx.values)) / minimum(abs.(F_Txx.values))

Tyy = assemble(T,Y,Y)
F_Tyy = eigen(Tyy)
@show maximum(abs.(F_Tyy.values)) / minimum(abs.(F_Tyy.values))

Nxy = assemble(N,X,Y)
@show cond(Nxy)

iNxy = inv(Nxy)
P = iNxy' * Tyy * iNxy
A = P * Txx

F_A = eigen(A)
@show maximum(abs.(F_A.values)) / minimum(abs.(F_A.values))
scatter(real.(F_A.values), imag.(F_A.values), markersize=0.01, limits=FRect(-1,-1,2,2))

E = Maxwell3D.planewave(direction=ẑ, polarization=ŷ, wavenumber=κ)
e = (n × E) × n
b = assemble(e,X)

solver = BEAST.GMRESSolver(A)
u, conv_hist = solve(solver,P*b)

Φ = range(0,stop=0,length=1)
Θ = range(0,stop=π,length=100)
pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]

F = Maxwell3D.farfield(wavenumber=κ)
ffd = potential(F, pts, u, X)
fcr, geo = facecurrents(u, X)

xs = range(-2,stop=2,length=51)
zs = range(-8,stop=4,length=150)
gridpoints = [point(x,0,z) for x in xs, z in zs]

N = MWSingleLayerField3D(wavenumber = κ)
nfd = potential(N, gridpoints, u, X)
nfd = reshape(nfd, size(gridpoints))
nfd .-= E.(gridpoints)

sc = Scene()
p1 = Makie.mesh(Γ,color=fc(norm.(fcr)),shading=false);
p2 = scatter(Θ, real.(norm.(ffd)),legend=false);
p3 = heatmap(clamp.(real.(norm.(nfd')), 0.0, 2.0));
p4 = lines(real.(getindex.(nfd,2))[25,:],legend=false);
s = vbox(hbox(p2,p1),hbox(p4,p3))

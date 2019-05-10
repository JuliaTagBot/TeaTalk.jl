using TeaTalk
d = pathof(TeaTalk)
Γ = readmesh(joinpath(dirname(d),"../examples/meshes","sphere2.in"))
X = raviartthomas(Γ)

κ = 1.0
sl = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=ŷ, wavenumber=κ)
e = (n × E) × n

@hilbertspace j
@hilbertspace k
efie = @discretise sl[k,j]==e[k]  j∈X k∈X
u = gmres(efie)

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

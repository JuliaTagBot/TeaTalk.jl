using TeaTalk

Γ = readmesh(joinpath(@__DIR__,"meshes","sphere2.in"))
X = raviartthomas(Γ)
Y = buffachristiansen(Γ)

κ = 1.0
μ = ϵ = 1.0
ω = κ/√(ϵ*μ)

dl = Maxwell3D.doublelayer(wavenumber=κ)
nx = BEAST.NCross()

E = Maxwell3D.planewave(direction=ẑ, polarization=ŷ, wavenumber=κ)
H = -1/(im*μ*ω)*curl(E)
h = (n × H) × n

@hilbertspace j
@hilbertspace k
mfie = @discretise (dl+0.5nx)[k,j] == h[k]  j∈X k∈Y
u = gmres(mfie)

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
p3 = heatmap(clamp.(real.(norm.(nfd')), 0.0, 2.0), colorbar=true);
p4 = lines(real.(getindex.(nfd,2))[25,:],legend=false);
s = vbox(hbox(p2,p1),hbox(p4,p3))

using TeaTalk
d = pathof(TeaTalk)
Γ = readmesh(joinpath(dirname(d),"../examples/meshes","sphere2.in"))
Y = buffachristiansen(Γ)

fcr1, geo1 = facecurrents(uv(123,numfunctions(Y)),Y)

scene = Makie.mesh(geo1, color=fc(norm.(fcr1)), shading=false, colormap=:RdYlBu)
wireframe!(scene[end][1], color = (:black, 0.6), linewidth = 3)
cam = cameracontrols(scene)
cam.eyeposition[] = Float32[1.67307, -1.11678, -0.659771]
cam.lookat[] = Float32[0.159456, 0.136933, 0.0]
cam.upvector[] = Float32[0.245085, -0.203002, 0.94801]
update_cam!(scene, cam)
scene.center = false
display(scene)

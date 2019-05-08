function meshdisk(;radius, h, udim=3)

    delta = h

    slope = 0.0
    dist = radius
    nRadii = ceil(Int, dist/delta)+1
    dStep = dist / (nRadii-1)
    rStep = dStep * radius / dist
    radii = range(0, step=rStep, stop=radius)

    circumf = 2*pi*radii
    nPhi    = ceil.(Int, circumf/delta)
    nPhi    = max.(nPhi,6)
    nPhi[1] = 1
    phiStep = 2*pi ./ nPhi
    offsets = [1; cumsum(nPhi).+1]

    closest_idx(phi,n) = mod1(round(Int,phi/2/pi*n)+1,n)
    gidx(i,j) = offsets[i] + mod1(j,nPhi[i])-1;

    NoV = sum(nPhi)
    V = zeros(NoV,3)
    V[1,:] = [0 0 0]
    for i = 1 : nRadii
        # phis = transpose(0 : phiStep(i) : 2*pi-phiStep(i));
        # phis = 0 : phiStep[i] : 2*pi-phiStep[i]
        phis = range(0, step=phiStep[i], length=nPhi[i])
        @assert length(phis) == nPhi[i]
        V[gidx(i,1) : gidx(i,1) + nPhi[i] - 1, :] = radii[i] * [cos.(phis) sin.(phis) slope*ones(nPhi[i],1)]
    end

    F = zeros(2*(NoV-1),3)
    idx = 1
    for i = 2 : nRadii
        phis = range(0, step=phiStep[i], length=nPhi[i])

        kprev = closest_idx(phis[end], nPhi[i-1])
        for j = 1 : nPhi[i]
            k = closest_idx(phis[j], nPhi[i-1]);

            if k != kprev
                F[idx,:] = [ gidx(i,j) gidx(i-1,k) gidx(i,j-1) ]; idx = idx + 1
                F[idx,:] = [ gidx(i-1,k) gidx(i-1,k-1) gidx(i,j-1)]; idx = idx + 1
            else
                F[idx,:] = [ gidx(i,j) gidx(i-1,k) gidx(i,j-1) ]; idx = idx + 1
            end
            kprev = k
        end
    end

    sV = [point(V[i,:]...) for i in axes(V,1)]
    sF = [index(F[i,:]...) for i in 1:idx-1]

    CompScienceMeshes.Mesh(sV,sF)
end

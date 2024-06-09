"""
    loglikelihoodargo(x, y, t, lonf, latf)

Calculate log-likelihood for Argo profiles with parameters x, observations y, 
times tf, longitudes lonf, and latitudes latf based on pointwise covariances. 
"""
function loglikelihoodargo(x, y, tf, lonf, latf; σtrend=0, σannual=0, itc0=0, tλx=200)

    # trend prior (mK/d)
    σtrend /= meanyear
    
    # real time (days)
    trd = Dates.value.(tf - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24

    # mid time
    tm = (min(trd...)+max(trd...))/2

    # parameters
    λt,λx,σs,σm = exp.(x)

    # trend prior scaling variance
    tc0 = σtrend^2/itc0
    
    # uniform seasonality covariance 
    Ξ = Diagonal((σannual*ones(4)/3.2).^2/2)

    # initialize log-likelihood gradient
    Δll = zeros(length(x))
    
    # time length
    m = length(trd)
    
    # seasonal frequency
    ω = 2π/SOT.meanyear
    
    # matrix to add seasonality
    A = [I(m) cos.(ω*trd) sin.(ω*trd) cos.(2ω*trd) sin.(2ω*trd)]
    
    # distance matrix
    dx = dist.(lonf', latf', lonf, latf)/1e3    

    # stochastic covariance 
    C0 = σs^2*sqrt(2)*exp.(-abs.(trd.-trd')/λt-dx/λx/sqrt(2))
    Cs = C0.*cos.(dx/λx/sqrt(2).-π/4)
    
    # offset and trend covariances
    Co,Ct = σm^2*exp.(-dx/tλx),(tc0*((trd.-tm)*(trd'.-tm))).*exp.(-dx/tλx)
    
    # full covariance without seasonality
    C = Cs+Co+Ct

    # add seasonality 
    Ryy = A*[C zeros(m, 4); zeros(4, m) Ξ]*A' 
    
    # invert covariance
    iRyy = inv(Matrix(Ryy))
    
    # log-likelihood
    ll = -0.5*(logdet(Ryy)+y'*iRyy*y+m*log(2π))

    # define gradient function
    α = iRyy*y
    δll0(δRyy) = 0.5*tr((α*α'-iRyy)*δRyy)

    # gradient along log-timescale
    δR1 = (abs.(trd.-trd')/λt).*Cs
    Δll[1] += δll0(δR1)
    
    # gradient along log-lengthscale
    δR2 = (dx/λx/sqrt(2)).*(Cs + C0.*sin.(dx/λx/sqrt(2).-π/4)) 
    Δll[2] += δll0(δR2)
    
    # gradient along log-stochastic scale
    δR3 = 2*(Cs) 
    Δll[3] += δll0(δR3)
    
    # gradient along log-offset scale
    δR4 = 2*(Co) 
    Δll[4] += δll0(δR4)

    return ll,Δll
end

"""
    cint(Δθ)

Calculate path-path stochastic covariance with a differential azimuth Δθ 
based on pointwise stochastic covariance. 
"""
function cint(Δθ;λx=65e-3,L=3.2,rtol=1e-6,atol=1e-6)

    # pointwise stochastic covariance 
    c0(x) = sqrt(2)*exp.(-sqrt(x[1]^2+x[2]^2-2*x[1]*x[2]*cosd(Δθ))/λx/sqrt(2)).*cos.(sqrt(x[1]^2+x[2]^2-2*x[1]*x[2]*cosd(Δθ))/λx/sqrt(2).-π/4)
    
    return hcubature(c0, [0,0], [L,L]; rtol=rtol, atol=atol)[1]
end

"""
    cinta(xp,yp,θp)

Calculate point-path stochastic covariance with point location (xp,yp) and path
azimuth θp based on pointwise stochastic covariance. 
"""
function cinta(xp,yp,θp;tstation=[166.90986,19.71786],event = [142.85,38.10],λx=65e-3,L=3.2,rtol=1e-6,atol=1e-6)

    # reference azimuth
    θ0 = azimuth(tstation[2],tstation[1],event[2],event[1])

    # differential azimuth
    Δθ = azimuth(tstation[2],tstation[1],yp,xp) - θ0 - θp

    # point-receiver distance    
    d = dist(xp, yp, tstation[1], tstation[2])/1e6
    
    # pointwise stochastic covariance
    c0(x) = sqrt(2)*exp.(-sqrt(x[1]^2+d^2-2x[1]*d*cosd(Δθ))/λx/sqrt(2)).*cos.(sqrt(x[1]^2+d^2-2x[1]*d*cosd(Δθ))/λx/sqrt(2).-π/4)
    
    return hcubature(c0, [0], [L]; rtol=rtol, atol=atol)[1]
end

"""
    cint2p(p1,p2,θ1,θ2)

Calculate stochastic covariance between two paths to two receiver locations 
p1, p2 with azimuths θ1, θ2 based on pointwise stochastic covariance. 
"""
function cint2p(p1,p2,θ1,θ2;event = [142.85,38.10],λx=65e-3,L=3.2,rtol=1e-6,atol=1e-6)
    
    # first azimuth to due north
    θ1 += azimuth(p1[2],p1[1],event[2],event[1])
    
    # second azimuth to due north
    θ2 += azimuth(p2[2],p2[1],event[2],event[1])
    
    # distance betwen receivers
    d = dist(p1[1], p1[2], p2[1], p2[2])/1e6
    
    # projected distances
    x20,y20 = d*cosd(θ2),d*sind(θ2)
    
    # distance in covariance integration
    cd(x) = sqrt((x[1]*cosd(θ1)-x[2]*cosd(θ2)-x20)^2+(x[1]*sind(θ1)-x[2]*sind(θ2)-y20)^2)
    
    # piontwise stochastic covariance
    c0(x) = sqrt(2)*exp.(-cd(x)/λx/sqrt(2)).*cos.(cd(x)/λx/sqrt(2).-π/4)
    
    return hcubature(c0, [0,0], [L,L]; rtol=rtol, atol=atol)[1]
end

"""
    cint2d(θp)

Calculate path-area stochastic covariance with a relative azimuth θp based on
pointwise stochastic covariance. 
"""
function cint2d(θp;shape="tri",θ1=-6.5,θ2=3.5,λx=65e-3,L=3.2,rtol=1e-6,atol=1e-6)   

    # calculation for a circular sector
    if shape=="tri"
    
        # integration limit
        θmax,θmin = [θ2,θ1] .- θp
        
        # pointwise stochastic covariance
        c0(x) = sqrt(2)*x[1]*exp(-sqrt(x[1]^2+x[2]^2-2*x[1]*x[2]*cosd(x[3]))/λx/sqrt(2))*cos(sqrt(x[1]^2+x[2]^2-2x[1]*x[2]*cosd(x[3]))/λx/sqrt(2)-π/4)
        
        return hcubature(c0, [0,0,θmin], [L,L,θmax]; rtol=rtol, atol=atol)[1]/((θmax-θmin)*L/2)
        
    # calculation for a rectangular sector
    elseif shape=="rec"
    
        # pointwise stochastic covariance
        c1(x) = sqrt(2)*exp.(-sqrt((x[1]*cosd(θp)-x[2])^2+(x[1]*sind(θp)-x[3])^2)/λx/sqrt(2)).*cos.(sqrt((x[1]*cosd(θp)-x[2])^2+(x[1]*sind(θp)-x[3])^2)/λx/sqrt(2).-π/4)
        
        return hcubature(c1, [0,0,π*L*θ1/180], [L,L,π*L*θ2/180]; rtol=rtol, atol=atol)[1]/(π*L*(θ2-θ1)/180)
    end
end

"""
    cint2da(xp,yp)

Calculate point-area stochastic covariance with point location (xp,yp) based on
pointwise stochastic covariance. 
"""
function cint2da(xp,yp;tstation=[166.90986,19.71786],event = [142.85,38.10],shape="tri",θ1=-6.5,θ2=3.5,λx=65e-3,L=3.2,rtol=1e-6,atol=1e-6) 

    # reference azimuth
    θ0 = azimuth(tstation[2],tstation[1],event[2],event[1])

    # differential azimuth
    Δθ = azimuth(tstation[2],tstation[1],yp,xp) - θ0 - θp

    # point-receiver distance    
    d = dist(xp, yp, tstation[1], tstation[2])/1e6
    
    
    # calculation for a circular sector
    if shape=="tri"
    
        # integration limit
        θmax,θmin = [θ2,θ1] .- θp
        
        # pointwise stochastic covariance
        c0(x) = sqrt(2)*x[1]*exp(-sqrt(x[1]^2+d^2-2*x[1]*d*cosd(x[2]))/λx/sqrt(2))*cos(sqrt(x[1]^2+d^2-2x[1]*d*cosd(x[2]))/λx/sqrt(2)-π/4)
        
        return hcubature(c0, [0,θmin], [L,θmax]; rtol=rtol, atol=atol)[1]/((θmax-θmin)*L/2)
        
    # calculation for a rectangular sector
    else
        
        # pointwise stochastic covariance
        c1(x) = sqrt(2)*exp.(-sqrt((d*cosd(θp)-x[1])^2+(d*sind(θp)-x[2])^2)/λx/sqrt(2)).*cos.(sqrt((d*cosd(θp)-x[1])^2+(d*sind(θp)-x[2])^2)/λx/sqrt(2).-π/4)
        
        return hcubature(c1, [0,π*L*θ1/180], [L,π*L*θ2/180]; rtol=rtol, atol=atol)[1]/(π*L*(θ2-θ1)/180)
    end
end

"""
    tcint(Δθ)

Calculate path-path trend covariance with a differential azimuth Δθ 
based on pointwise trend covariance. 
"""
function tcint(Δθ;λx=0.2,L=3.2,rtol=1e-6,atol=1e-6) 

    # pointwise trend covariance
    c0(x) = exp.(-sqrt(x[1]^2+x[2]^2-2x[1]*x[2]*cosd(Δθ))/λx)
    
    return hcubature(c0, [0,0], [L,L]; rtol=rtol, atol=atol)[1]
end

"""
    tcinta(xp,yp,θp)

Calculate point-path trend covariance with point location (xp,yp) and path
azimuth θp based on pointwise trend covariance. 
"""
function tcinta(xp,yp,θp;tstation=[166.90986,19.71786],event = [142.85,38.10],λx=0.2,L=3.2,rtol=1e-6,atol=1e-6) 
    # reference azimuth
    θ0 = azimuth(tstation[2],tstation[1],event[2],event[1])

    # differential azimuth
    Δθ = azimuth(tstation[2],tstation[1],yp,xp) - θ0 - θp

    # point-receiver distance    
    d = dist(xp, yp, tstation[1], tstation[2])/1e6
    
    # pointwise trend covariance
    c0(x) = exp.(-sqrt(x[1]^2+d^2-2x[1]*d*cosd(Δθ))/λx)
    
    return hcubature(c0, [0], [L]; rtol=rtol, atol=atol)[1]
end

"""
    tcint2d(θp)

Calculate path-area trend covariance with a relative azimuth θp based on
pointwise trend covariance. 
"""
function tcint2d(θp;shape="tri",θ1=-6.5,θ2=3.5,λx=0.2,L=3.2,rtol=1e-6,atol=1e-6) 

    # calculation for a circular sector
    if shape=="tri"
    
        # integration limit
        θmax,θmin = [θ2,θ1] .- θp
        
        # pointwise trend covariance
        c0(x) = x[1]*exp(-sqrt(x[1]^2+x[2]^2-2*x[1]*x[2]*cosd(x[3]))/λx)
        
        return hcubature(c0, [0,0,θmin], [L,L,θmax]; rtol=rtol, atol=atol)[1]/((θmax-θmin)*L/2)
        
    # calculation for a triangular sector
    elseif shape=="rec"
    
        # pointwise trend covariance
        c1(x) = exp.(-sqrt((x[1]*cosd(θp)-x[2])^2+(x[1]*sind(θp)-x[3])^2)/λx)
        
        return hcubature(c1, [0,0,π*L*θ1/180], [L,L,π*L*θ2/180]; rtol=rtol, atol=atol)[1]/(π*(θ2-θ1)*L/180)
    end
end

"""
    tcint2da(xp,yp)

Calculate point-area trend covariance with point location (xp,yp) based on
pointwise trend covariance. 
"""
function tcint2da(xp,yp;tstation=[166.90986,19.71786],event = [142.85,38.10],shape="tri",θ1=-6.5,θ2=3.5,λx=0.2,L=3.2,rtol=1e-6,atol=1e-6) 

    # reference azimuth
    θ0 = azimuth(tstation[2],tstation[1],event[2],event[1])

    # differential azimuth
    Δθ = azimuth(tstation[2],tstation[1],yp,xp) - θ0 - θp

    # point-receiver distance    
    d = dist(xp, yp, tstation[1], tstation[2])/1e6
    
    # calculation for a circular sector
    if shape=="tri"
    
        # integration limit
        θmax,θmin = [θ2,θ1] .- θp
        
        # pointwise trend covariance
        c0(x) = x[1]*exp(-sqrt(x[1]^2+d^2-2*x[1]*d*cosd(x[2]))/λx)
        
        return hcubature(c0, [0,θmin], [L,θmax]; rtol=rtol, atol=atol)[1]/((θmax-θmin)*L/2)
    
    # calculation for a triangular sector
    else
    
        # pointwise trend covariance
        c1(x) = exp.(-sqrt((d*cosd(θp)-x[1])^2+(d*sind(θp)-x[2])^2)/λx)
        
        return hcubature(c1, [0,π*L*θ1/180], [L,π*L*θ2/180]; rtol=rtol, atol=atol)[1]/(π*(θ2-θ1)*L/180)
    end
end

"""
find nearby argo profiles
"""
function nearargo(floats,events,nxk,d0,Δt,tstation)
    m = size(events,1)
    nearfloats = []
    for i = 1:m

        datestart = events[i,:time] - Dates.Day(Δt)
        dateend = events[i,:time] + Dates.Day(Δt)
        xi = events[i,:longitude]#0.5*(ppairs[i,:longitude1]+ppairs[i,:longitude2]) 
        yi = events[i,:latitude]#0.5*(ppairs[i,:latitude1]+ppairs[i,:latitude2]) 
    
        xki, yki = SOT.findpath(tstation, (xi,yi), nxk)
        dateidx = (datestart .≤ floats.t .≤ dateend) #.|| (datestart2 .≤ floats.t .≤ dateend2)
        floatsi = floats[dateidx,:]
    
        distancei = minimum(SOT.dist.(xki', yki', floatsi.x, floatsi.y),dims=2)./1e3
        
        if i == 1
            nearfloats = floatsi[vec(distancei) .≤ d0,:]
        else
            nearfloats = vcat(nearfloats, floatsi[vec(distancei) .≤ d0,:], cols = :union)
            unique!(nearfloats, [:x,:y,:t])
        end
    
        @printf("%d events use %d profiles\n",i,size(nearfloats,1))
    end
    return nearfloats
end
"""
linvert travel time anomalies using Argo and T-wave data
"""
function invertargot0(te,azme,x,σtrend, σannual, σsemiannual,tpairs,ppairs,tstation,evtpos,pstations,nearfloats,itp_c,itp_ca,itp_tc,itp_tca;tλx=200)
    ### twave
    trde = Dates.value.(te - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
    tm = (max(trde...)+min(trde...))/2
    m = length(te)
    
    ### twave
    _,_,θd = getazimuth(tpairs,ppairs,tstation[1],tstation[2],evtpos[1],evtpos[2])
    t1,E1,Ecs,Esn,Ecsh,Esnh = SOT.getE(tpairs,ppairs,tstation,pstations;hydro=true)

    nt,np = size(tpairs,1),size(ppairs,1)
    σtrend /= meanyear
    tc0 = σtrend^2/itp_tc(0)
    λt,λx,σs,σm,σn,σp,σx,σ,σnp,σh = x
    
    l = 1
    m1 = length(t1)
    # real time (days)
    trd = Dates.value.(t1 - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
    C = σs^2*exp.(-abs.(trd.-trd')/λt).*itp_c.(abs.((θd.-θd')))
    C += tc0*((trd.-tm)*(trd'.-tm)).*itp_tc.(abs.((θd.-θd')))
    E1 = [E1[:,1:2m1] E1[:,2m1+2:end]]

    # covariance matrix assuming no correlation between singular vector expansion coefficients
    R = [C zeros(l*m1, m1); zeros(m1, (l+1)*m1)] + σp^2*kron(sparse(ones(l+1, l+1)), I(m1))
    #R = [R zeros(size(R, 1), l); zeros(l, size(R, 2)) (σtrend.^2)]
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) I(2)*(σannual.^2/2)]
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) I(2)*(σsemiannual.^2/2)]

    # noise covariance
    Rcsn = Ecs*Ecs'+Esn*Esn'
    Rcsnh = Ecsh*Ecsh'+Esnh*Esnh'
    Nt = σ^2*spdiagm([ones(nt);zeros(np)])+σx^2*Rcsn+σh^2*Rcsnh
    Nt += σnp^2*spdiagm([zeros(nt);ones(np)])

    Ct = E1*R*E1'+Nt

    N = σn^2*I

    C0 = σs^2*exp.(-abs.(trde.-trde')/λt).*itp_c.(abs.((azme.-azme')))

    C0 = [C0 zeros(m, m+4); zeros(m, m) tc0*itp_tc.(abs.((azme.-azme'))) zeros(m,4);zeros(4,2m) R[end-3:end,end-3:end]]
    
    Rt = σs^2*exp.(-abs.(trd.-trde')/λt).*itp_c.(abs.(θd.-azme')) 

    Rt = [Rt tc0*(trd.-tm).*itp_tc.(abs.(θd.-azme')) zeros(m1,4)]
    Rt = [Rt; zeros(m1,2m+4); zeros(4,2m) R[end-3:end,end-3:end]]
    Rt = E1*Rt
    
    @printf("path covariance finished, time: %s\n",Dates.format(now(), "HH:MM:SS"))

    trf = Dates.value.(nearfloats.t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
    ma = length(trf)
    dx = SOT.dist.(nearfloats.x', nearfloats.y', nearfloats.x, nearfloats.y)/1e3
    Ca = σs^2*sqrt(2)*exp.(-abs.(trf.-trf')/λt-dx/λx/sqrt(2))
    Ca = Ca.*cos.(dx/λx/sqrt(2).-π/4)
    Ca += tc0*((trf.-tm)*(trf'.-tm)).*exp.(-dx/tλx)
    Ca += σm^2*exp.(-dx/tλx)
    Ca = [Ca zeros(ma, 4); zeros(4, ma) R[end-3:end,end-3:end]/3.2^2]

    Ra = σs^2*exp.(-abs.(trf.-trde')/λt).*itp_ca.(nearfloats.x, nearfloats.y, azme')
    Ra = [Ra tc0*(trf.-tm).*itp_tca.(nearfloats.x, nearfloats.y, azme')]
    ω = 2π/SOT.meanyear
    Ra = [Ra σannual.^2/2*cos.(ω*trf)/3.2 σannual.^2/2*sin.(ω*trf)/3.2] 
    Ra = [Ra σsemiannual.^2/2*cos.(2ω*trf)/3.2 σsemiannual.^2/2*sin.(2ω*trf)/3.2]

    Ea = [I(ma) cos.(ω*trf) sin.(ω*trf) cos.(2ω*trf) sin.(2ω*trf)]
    Ca = Ea*Ca*Ea'+N

    @printf("point covariance finished, time: %s\n",Dates.format(now(), "HH:MM:SS"))
    
    iCa = inv(Array(Ca))
    Pa = C0-Ra'*iCa*Ra
    xa = Ra'*iCa*nearfloats.Δsp1

    @printf("point inversion finished, time: %s\n",Dates.format(now(), "HH:MM:SS"))

    iCt = inv(Array(Ct))
    Pt = C0-Rt'*iCt*Rt
    xt = Rt'*iCt*[tpairs.Δτ; ppairs.Δτ]
    
    @printf("path covariance finished, time: %s\n",Dates.format(now(), "HH:MM:SS"))

    Rat = σs^2*exp.(-abs.(trf.-trd')/λt).*itp_ca.(nearfloats.x, nearfloats.y, θd')
    Rat += tc0*((trf.-tm)*(trd'.-tm)).*itp_tca.(nearfloats.x, nearfloats.y, θd')
    Rat = [Rat zeros(ma,m1) Ra[:,end-3:end]]
    Rat = Rat*E1'
    
    @printf("combined covariance finished, time: %s\n",Dates.format(now(), "HH:MM:SS"))

    RiPat = (Ra'-Rt'*iCt*Rat')*inv(Ca-Rat*iCt*Rat')
    RiPta = (Rt'-Ra'*iCa*Rat)*inv(Ct-Rat'*iCa*Rat)
    xat = RiPat*nearfloats.Δsp1 + RiPta*[tpairs.Δτ; ppairs.Δτ]
    Pat = C0-(RiPat*Ra+RiPta*Rt)
    
    @printf("combined inversion finished, time: %s\n",Dates.format(now(), "HH:MM:SS"))
    
    dC = Ra'*iCa*Rat*iCt*Rt 

    return Pa,xa,Pt,xt,Pat,xat,C0,dC
end

"""
linvert ship ctd using Argo
"""
function invertctd(σtrend, σannual,ctds,nearfloats,λx,λt,σs,σn,itp_tc;tλx=200)
    trc = Dates.value.(ctds.t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
    tm = (max(trc...)+min(trc...))/2
    m = length(trc)

    N = σn^2*I
    σtrend /= meanyear
    tc0 = σtrend^2/itp_tc(0)
    Ξ = Diagonal((σannual*ones(4)/sqrt(2)/3.2).^2)
    
    dxc = SOT.dist.(ctds.x', ctds.y', ctds.x, ctds.y)/1e3
    C0 = σs^2*sqrt(2)*exp.(-abs.(trc.-trc')/λt-dxc/λx/sqrt(2))
    C0 = C0.*cos.(dxc/λx/sqrt(2).-π/4)
    C0 += tc0*((trc.-tm)*(trc'.-tm)).*exp.(-dxc/tλx)

    C0 = [C0 zeros(m, 4); zeros(4,m) Ξ]

    trf = Dates.value.(nearfloats.t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
    ma = length(trf)
    dx = SOT.dist.(nearfloats.x', nearfloats.y', nearfloats.x, nearfloats.y)/1e3
    Ca = σs^2*sqrt(2)*exp.(-abs.(trf.-trf')/λt-dx/λx/sqrt(2))
    Ca = Ca.*cos.(dx/λx/sqrt(2).-π/4)
    Ca += tc0*((trf.-tm)*(trf'.-tm)).*exp.(-dx/tλx)
    Ca = [Ca zeros(ma, 4); zeros(4, ma) Ξ]

    dxa = SOT.dist.(ctds.x', ctds.y', nearfloats.x, nearfloats.y)/1e3
    Ra = σs^2*sqrt(2)*exp.(-abs.(trf.-trc')/λt-dxa/λx/sqrt(2))
    Ra = Ra.*cos.(dxa/λx/sqrt(2).-π/4)
    Ra += tc0*((trf.-tm)*(trc'.-tm)).*exp.(-dxa/tλx)
    Ra = [Ra zeros(ma,4); zeros(4,m) Ξ]
    ω = 2π/SOT.meanyear

    E = [I(m) cos.(ω*trc) sin.(ω*trc) cos.(2ω*trc) sin.(2ω*trc)]
    C0 = E*C0*E'+N
    Ea = [I(ma) cos.(ω*trf) sin.(ω*trf) cos.(2ω*trf) sin.(2ω*trf)]
    Ca = Ea*Ca*Ea'+N
    Ra = Ea*Ra*E'

    iCa = inv(Matrix(Ca))
    Pa = C0-Ra'*iCa*Ra
    xa = Ra'*iCa*nearfloats.Δsp1

    return Pa,xa
end

"""
linvert travel time anomalies using Argo and T-wave data
"""
function invertargottc(te,azme,x,σtrend, σannual, σsemiannual,tpairs,ppairs,tstation,evtpos,pstations,nearfloats;
                       itp_c=cint,itp_ca=cinta,itp_tc=tcint,itp_tca=tcinta,tλx=200,itp=true,mode="self")#tσs=0.31,tλt=500,
    #vtrend=true
    λt,λx,σs,σm,σn,σp,σx,σ,σnp,σh = x
    println(itp_c)
    
    trde = Dates.value.(te - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
    tm = (max(trde...)+min(trde...))/2
    m = length(te)
    
    σtrend /= meanyear
    tc0 = σtrend^2/tcint(0;λx=λx/1e3)
    
    N = σn^2*I
    l = 1

    if itp
        C0 = itp_c.(abs.(azme.-azme'))
        Ct0 = itp_tc.(abs.(azme.-azme'))
    else
        C0 = ones(m,m)*itp_c(0;λx=λx/1e3)
        Ct0 = ones(m,m)*itp_tc(0;λx=λx/1e3)
        if mode!="ones"
            for i = 1:m
              for j = i+1:m
                C0[i,j] = itp_c(abs(azme[i]-azme[j]);λx=λx/1e3)
                Ct0[i,j] = itp_tc(abs(azme[i]-azme[j]);λx=λx/1e3)
                C0[j,i],Ct0[j,i] = C0[i,j],Ct0[i,j]
              end
            end
        end
    end
    C0 = σs^2*exp.(-abs.(trde.-trde')/λt).*C0
    
    nfit = 3
    ftc(τ1,τ2) = tc0*((τ1 .- tm)*(τ2 .- tm)) #: tσs^2*sqrt(2)*exp.(-abs.(τ1 .- τ2)/tλt/sqrt(2)).*cos.(abs.(τ1 .- τ2)/tλt/sqrt(2).-π/4)

    C0 = [C0 zeros(m,m); zeros(m,m) tc0*Ct0]

    ### twave
    _,_,θd = getazimuth(tpairs,ppairs,tstation[1],tstation[2],evtpos[1],evtpos[2])
    t1,E1,Ecs,Esn,Ecsh,Esnh = SOT.getE(tpairs,ppairs,tstation,pstations;hydro=true)
    
    nt,np = size(tpairs,1),size(ppairs,1)
    m1 = length(t1)
    # real time (days)
    trd = Dates.value.(t1 - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
    if itp
        C = σs^2*exp.(-abs.(trd.-trd')/λt).*itp_c.(abs.(θd.-θd'))
        C += ftc(trd,trd').*itp_tc.(abs.(θd.-θd'))
    else
        if mode == "self"
            C = C0[1:m,1:m]+ftc(trd,trd').*C0[m+1:2m,m+1:2m]/tc0
        else
            C = ones(m1,m1)*itp_c(0;λx=λx/1e3)
            Ct = ones(m1,m1)*itp_tc(0;λx=λx/1e3)
            for i = 1:m1
              for j = i+1:m1
                C[i,j] = itp_c(abs(θd[i]-θd[j]);λx=λx/1e3)
                Ct[i,j] = itp_tc(abs(θd[i]-θd[j]);λx=λx/1e3)
                C[j,i],Ct[j,i] = C[i,j],Ct[i,j]
              end
            end
            C = σs^2*exp.(-abs.(trd.-trd')/λt).*C
            C += ftc(trd,trd').*Ct
        end
    end

    # covariance matrix assuming no correlation between singular vector expansion coefficients
    R = [C zeros(l*m1, m1); zeros(m1, (l+1)*m1)] + σp^2*kron(sparse(ones(l+1, l+1)), I(m1))
    E1 = [E1[:,1:2m1] E1[:,2m1+2:end]]

    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) I(2)*(σannual.^2/2)]
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) I(2)*(σsemiannual.^2/2)]

    # noise covariance
    Rcsn = Ecs*Ecs'+Esn*Esn'
    Rcsnh = Ecsh*Ecsh'+Esnh*Esnh'
    Nt = σ^2*spdiagm([ones(nt);zeros(np)])+σx^2*Rcsn+σh^2*Rcsnh
    Nt += σnp^2*spdiagm([zeros(nt);ones(np)])

    Ct = E1*R*E1'+Nt
    
    if itp
        Rt = σs^2*exp.(-abs.(trd.-trde')/λt).*itp_c.(abs.(θd.-azme')) 
        Rt = [Rt tc0*(trd.-tm).*itp_tc.(abs.(θd.-azme'))]
    else
        if mode=="self"
            Rt = [C0[1:m,1:m] (trd.-tm).*C0[m+1:2m,m+1:2m]]
        elseif mode=="ones"
            Rt = σs^2*exp.(-abs.(trd.-trde')/λt).*repeat(itp_c.(abs.((θd.-azme[1]));λx=λx/1e3),1,m)
            Rt = [Rt tc0*(trd.-tm).*repeat(itp_tc.(abs.(θd.-azme[1]);λx=λx/1e3),1,m)]
        else
            Rt = σs^2*exp.(-abs.(trd.-trde')/λt).*itp_c.(abs.((θd.-azme'));λx=λx/1e3) 
            Rt = [Rt tc0*(trd.-tm).*itp_tc.(abs.(θd.-azme');λx=λx/1e3)]
        end
    end
    Rt = [Rt zeros(m1,nfit+1); zeros(m1,size(C0,1)+nfit+1); zeros(nfit+1,size(C0,2)) R[end-nfit:end,end-nfit:end]]
    Rt = E1*Rt
    @printf("T wave covariance finished, time: %s\n",Dates.format(now(), "HH:MM:SS"))
    
    C0 = [C0 zeros(size(C0,1), nfit+1); zeros(nfit+1,size(C0,2)) R[end-nfit:end,end-nfit:end]]

    trf = Dates.value.(nearfloats.t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
    ma = length(trf)
    dx = SOT.dist.(nearfloats.x', nearfloats.y', nearfloats.x, nearfloats.y)/1e3
    Ca = σs^2*sqrt(2)*exp.(-abs.(trf.-trf')/λt-dx/λx/sqrt(2))
    Ca = Ca.*cos.(dx/λx/sqrt(2).-π/4)
    Ca += tc0*((trf.-tm)*(trf'.-tm)).*exp.(-dx/tλx)

    @printf("add mean covariance in point data.\n")
    #σm = σs
    Ca += σm^2*exp.(-dx/tλx)

    Ca = [Ca zeros(ma, nfit+1); zeros(nfit+1, ma) R[end-nfit:end,end-nfit:end]/3.2^2]
   
    if itp
        Ra = σs^2*exp.(-abs.(trf.-trde')/λt).*itp_ca.(nearfloats.x, nearfloats.y, azme')
        Ra = [Ra tc0*(trf.-tm).*itp_tca.(nearfloats.x, nearfloats.y, azme')]
    else
        if mode=="ones"
            Ra = σs^2*exp.(-abs.(trf.-trde')/λt).*repeat(itp_ca.(nearfloats.x, nearfloats.y, azme[1];λx=λx/1e3),1,m)
            Ra = [Ra tc0*(trf.-tm).*repeat(itp_tca.(nearfloats.x, nearfloats.y, azme[1];λx=λx/1e3),1,m)]
        else
            Ra = σs^2*exp.(-abs.(trf.-trde')/λt).*itp_ca.(nearfloats.x, nearfloats.y, azme';λx=λx/1e3)
            Ra = [Ra tc0*(trf.-tm).*itp_tca.(nearfloats.x, nearfloats.y, azme';λx=λx/1e3)]
        end
    end
    ω = 2π/SOT.meanyear
    Ra = [Ra σannual.^2/2*cos.(ω*trf)/3.2 σannual.^2/2*sin.(ω*trf)/3.2]
    Ea = [I(ma) cos.(ω*trf) sin.(ω*trf) cos.(2ω*trf) sin.(2ω*trf)]

    Ra = [Ra σsemiannual.^2/2*cos.(2ω*trf)/3.2 σsemiannual.^2/2*sin.(2ω*trf)/3.2]

    Ca = Ea*Ca*Ea'+N
    @printf("point data covariance finished, time: %s\n",Dates.format(now(), "HH:MM:SS"))

    iCa = inv(Matrix(Ca))
    Pa = C0-Ra'*iCa*Ra
    xa = Ra'*iCa*nearfloats.Δsp1
    @printf("point data inversion finished, time: %s\n",Dates.format(now(), "HH:MM:SS"))

    iCt = inv(Matrix(Ct))
    Pt = C0-Rt'*iCt*Rt
    xt = Rt'*iCt*[tpairs.Δτ; ppairs.Δτ]
    @printf("T wave inversion finished, time: %s\n",Dates.format(now(), "HH:MM:SS"))

    if itp
        Rat = σs^2*exp.(-abs.(trf.-trd')/λt).*itp_ca.(nearfloats.x, nearfloats.y, θd')
        Rat += ftc(trf,trd').*itp_tca.(nearfloats.x, nearfloats.y, θd')
    else
        if mode == "self"
            Rat = Ra[:,1:m] + (trd'-tm).*Ra[:,m+1:2m]
        else
            Rat = σs^2*exp.(-abs.(trf.-trd')/λt).*itp_ca.(nearfloats.x, nearfloats.y, θd';λx=λx/1e3)
            Rat += ftc(trf,trd').*itp_tca.(nearfloats.x, nearfloats.y, θd';λx=λx/1e3)
        end
    end
    Rat = [Rat zeros(ma,m1) Ra[:,end-nfit:end]]
    Rat = Rat*E1'

    RiPat = (Ra'-Rt'*iCt*Rat')*inv(Ca-Rat*iCt*Rat')
    RiPta = (Rt'-Ra'*iCa*Rat)*inv(Ct-Rat'*iCa*Rat)
    xat = RiPat*nearfloats.Δsp1 + RiPta*[tpairs.Δτ; ppairs.Δτ]
    Pat = C0-(RiPat*Ra+RiPta*Rt)
    @printf("joint inversion finished, time: %s\n",Dates.format(now(), "HH:MM:SS"))
    
    dC = Ra'*iCa*Rat*iCt*Rt 

    #iC = inv([Ca Rat;Rat' Ct])
    #Rat = [Ra;Rt]
    #Pat = C0-Rat'*iC*Rat
    #xat = Rat'*iC*[nearfloats.Δsp1; tpairs.Δτ; ppairs.Δτ]]

    #ΔP = 2*Rat'*iC*Rat-Ra'*iCa*Ra-Rt'*iCt*Rt
    return Pa,xa,Pt,xt,Pat,xat,C0,dC
end

"""
invert sub-basin average travel time anomalies using Argo and T waves
"""
function invertargot2d(te,x,σtrend, σannual, σsemiannual,tpairs,ppairs,tstation,evtpos,pstations,nearfloats,λx,σs,σn,nxk,interp_linear,itpat,itpt,itp,c0;full=false)
    ### twave
    _,_,θd = getazimuth(tpairs,ppairs,tstation[1],tstation[2],evtpos[1],evtpos[2])
    t1,E1,Ecs,Esn,Ecsh,Esnh = SOT.getE(tpairs,ppairs,tstation,pstations;hydro=true)

    nt,np = size(tpairs,1),size(ppairs,1)
    σtrend /= meanyear
    λt,λx,σs,σm,σn,σp,σx,σ,σnp,σh = x
    l = 1
    m1 = length(t1)
    # real time (days)
    trd = Dates.value.(t1 - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
    C = σs^2*exp.(-abs.(trd.-trd')/λt).*interp_linear.(abs.((θd.-θd')))

    # covariance matrix assuming no correlation between singular vector expansion coefficients
    R = [C zeros(l*m1, m1); zeros(m1, (l+1)*m1)] + σp^2*kron(sparse(ones(l+1, l+1)), I(m1))
    R = [R zeros(size(R, 1), l); zeros(l, size(R, 2)) (σtrend.^2)]
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) I(2)*(σannual.^2/2)]
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) I(2)*(σsemiannual.^2/2)]

    # noise covariance
    Rcsn = Ecs*Ecs'+Esn*Esn'
    Rcsnh = Ecsh*Ecsh'+Esnh*Esnh'
    Nt = σ^2*spdiagm([ones(nt);zeros(np)])+σx^2*Rcsn+σh^2*Rcsnh
    Nt += σnp^2*spdiagm([zeros(nt);ones(np)])

    Ct = E1*R*E1'+Nt

    N = 0*I

    trde = Dates.value.(te - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
    if full
        m = length(te)
        tm = (max(trde...)+min(trde...))/2

        C0 = σs^2*c0*exp.(-abs.(trde.-trde')/λt)
        C0 = [C0 zeros(m, 5); zeros(5, m) R[end-4:end,end-4:end]]

        trf = Dates.value.(nearfloats.t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
        ma = length(trf)
        dx = SOT.dist.(nearfloats.x', nearfloats.y', nearfloats.x, nearfloats.y)/1e3
        Ca = σs^2*sqrt(2)*exp.(-abs.(trf.-trf')/λt-dx/λx/sqrt(2))
        Ca = Ca.*cos.(dx/λx/sqrt(2).-π/4)
        Ca = [Ca zeros(ma, 5); zeros(5, ma) R[end-4:end,end-4:end]/3.2^2]

        Ra = σs^2*exp.(-abs.(trf.-trde')/λt).*itp.(nearfloats.x, nearfloats.y)
        Ra = [Ra σtrend.^2*(trf.-tm)/3.2]
        ω = 2π/SOT.meanyear
        Ra = [Ra σannual.^2/2*cos.(ω*trf)/3.2 σannual.^2/2*sin.(ω*trf)/3.2] 
        Ra = [Ra σsemiannual.^2/2*cos.(2ω*trf)/3.2 σsemiannual.^2/2*sin.(2ω*trf)/3.2]

        Ea = [I(ma) (trf.-tm) cos.(ω*trf) sin.(ω*trf) cos.(2ω*trf) sin.(2ω*trf)]
        Ca = Ea*Ca*Ea'+N

        iCa = inv(Matrix(Ca))
        Pa = C0-Ra'*iCa*Ra
        xa = Ra'*iCa*nearfloats.Δsp1
    
        Rt = σs^2*exp.(-abs.(trd.-trde')/λt).*itpt.(θd)

        Rt = [Rt zeros(m1,5)]
        Rt = [Rt; zeros(m1,m+5); zeros(5,m) R[end-4:end,end-4:end]]
        Rt = E1*Rt

        iCt = inv(Matrix(Ct))
        Pt = C0-Rt'*iCt*Rt
        xt = Rt'*iCt*[tpairs.Δτ; ppairs.Δτ]

        Rat = σs^2*exp.(-abs.(trf.-trd')/λt).*itpat.(nearfloats.x, nearfloats.y, θd')
        Rat = [Rat zeros(ma,m1) Ra[:,end-4:end]]
        Rat = Rat*E1'

        iC = inv(Matrix([Ca Rat;Rat' Ct]))
        Rat = [Ra;Rt]
        Pat = C0-Rat'*iC*Rat
        xat = Rat'*iC*[nearfloats.Δsp1; tpairs.Δτ; ppairs.Δτ]

        return Pa,xa,Pt,xt,Pat,xat

    else
        year0,year1 = max(2008,Dates.year(min(nearfloats.t...))),min(2021,Dates.year(max(nearfloats.t...)))

        iC,iCm = [],[]

        for year = year0:year1
            ### target
            idxe = DateTime(year, 1, 1) .≤ te .< DateTime(year+1, 1, 1)
            m = length(te[idxe])

            C0 = σs^2*c0*exp.(-abs.(trde[idxe].-trde[idxe]')/λt)
            iC0 = inv(C0)

            iR = [iC0 zeros(m, m1+5); zeros(m1+5, m) inv(R1[m1+1:end,m1+1:end])]
        
            Ct = σs^2*exp.(-abs.(trd.-trde[idxe]')/λt).*itpt.(θd)

            A = [Ct σp^2*I(m1) zeros(m1,5);zeros(m1,m) σp^2*I(m1) zeros(m1,5);zeros(5,m+m1) R1[end-4:end,end-4:end]]
        
            ### Argo
            idx = DateTime(year, 1, 1) .≤ nearfloats.t .< DateTime(year+1, 1, 1)

            floatsi = nearfloats[idx,:]
            trfloats = Dates.value.(floatsi.t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
            dxi = SOT.dist.(floatsi.x', floatsi.y', floatsi.x, floatsi.y)/1e3
            Ci = σs^2*sqrt(2)*exp.(-abs.(trfloats.-trfloats')/λt-dxi/λx/sqrt(2))
            Ci = Ci.*cos.(dxi/λx/sqrt(2).-π/4)
            mfi = size(floatsi,1)

            @printf("year %d has %d profiles\n",year,mfi)
            Ri = σs^2*exp.(-abs.(trfloats.-trde[idxe]')/λt).*itp.(floatsi.x, floatsi.y)
            Rti = σs^2*exp.(-abs.(trfloats.-trd')/λt).*itpat.(floatsi.x, floatsi.y, θd')

            Ai = [E1*A; Ri zeros(mfi, m1+5)]
            Σti = [E1[1:nt,1:m1]*Rti'; zeros(np,mfi)]
            Σti = [Rt Σti; Rti*E1[1:nt,1:m1]' zeros(mfi,np) Ci+N]-Ai*iR*Ai'

            T = [zeros(mfi,np+nt) I(mfi); I(nt+np) zeros(nt+np,mfi)]
            Σti = T*Σti*T'
            iΓi = inv(Array(Σti))

            push!(iC,iR+iR*Ai'*T'*iΓi*T*Ai*iR)
            push!(iCm,iR*Ai'*T'*iΓi*[floatsi.Δsp1; tpairs.Δτ; ppairs.Δτ])
        end

        return iC,iCm
    end
end

"""
invert sub-basin average travel time anomalies using Argo and T waves
"""
function invertargot2dtc(te,x,σtrend, σannual, σsemiannual,tpairs,ppairs,tstation,evtpos,pstations,nearfloats,c0,c02d;
                         itp_c=cint,itp_ca=cinta,itp_tc=tcint,itp_tca=tcinta,itp_c2d=cint2d,itp_ca2d=cint2da,itp_tc2d=tcint2d,itp_tca2d=tcint2da,tλx=200,full=false)
    trde = Dates.value.(te - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
    tm = (max(trde...)+min(trde...))/2

    ### twave
    _,_,θd = getazimuth(tpairs,ppairs,tstation[1],tstation[2],evtpos[1],evtpos[2])
    t1,E1,Ecs,Esn,Ecsh,Esnh = SOT.getE(tpairs,ppairs,tstation,pstations;hydro=true)

    nt,np = size(tpairs,1),size(ppairs,1)
    σtrend /= meanyear
    tc0 = σtrend^2/itp_tc(0)
    λt,λx,σs,σm,σn,σp,σx,σ,σnp,σh = x
    l = 1
    m1 = length(t1)
    # real time (days)
    trd = Dates.value.(t1 - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
    C = abs.(θd.-θd')
    for i = 1:size(C,1)
      for j = i:size(C,2)
        C[i,j] = σs^2*exp(-abs(trd[i]-trd[j])/λt)*itp_c(C[i,j]) + tc0*((trd[i]-tm)*(trd[j]-tm))*itp_tc(C[i,j])
        if j>i
          C[j,i] = C[i,j]
        end
      end
    end
       
    #C = σs^2*exp.(-abs.(trd.-trd')/λt).*itp_c.(abs.((θd.-θd')))
    #C += tc0*((trd.-tm)*(trd'.-tm)).*itp_tc.(abs.((θd.-θd')))
    E1 = [E1[:,1:2m1] E1[:,2m1+2:end]]

    # covariance matrix assuming no correlation between singular vector expansion coefficients
    R = [C zeros(l*m1, m1); zeros(m1, (l+1)*m1)] + σp^2*kron(sparse(ones(l+1, l+1)), I(m1))
    #R = [R zeros(size(R, 1), l); zeros(l, size(R, 2)) (σtrend.^2)]
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) I(2)*(σannual.^2/2)]
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) I(2)*(σsemiannual.^2/2)]

    # noise covariance
    Rcsn = Ecs*Ecs'+Esn*Esn'
    Rcsnh = Ecsh*Ecsh'+Esnh*Esnh'
    Nt = σ^2*spdiagm([ones(nt);zeros(np)])+σx^2*Rcsn+σh^2*Rcsnh
    Nt += σnp^2*spdiagm([zeros(nt);ones(np)])

    Ct = E1*R*E1'+Nt

    N = 0*I

    m = length(te)

    C0 = σs^2*c0*exp.(-abs.(trde.-trde')/λt)

    C0 = [C0 zeros(m, 5); zeros(1, m) c02d*tc0 zeros(1,4); zeros(4,1+m) R[end-3:end,end-3:end]]

    trf = Dates.value.(nearfloats.t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
    ma = length(trf)
    dx = SOT.dist.(nearfloats.x', nearfloats.y', nearfloats.x, nearfloats.y)/1e3
    Ca = σs^2*sqrt(2)*exp.(-abs.(trf.-trf')/λt-dx/λx/sqrt(2))
    Ca = Ca.*cos.(dx/λx/sqrt(2).-π/4)
    Ca += (σm^2 .+ tc0*((trf.-tm)*(trf'.-tm))).*exp.(-dx/tλx)
    Ca = [Ca zeros(ma, 4); zeros(4, ma) R[end-3:end,end-3:end]/3.2^2]

    Ra = σs^2*exp.(-abs.(trf.-trde')/λt).*itp_ca2d.(nearfloats.x, nearfloats.y)
    Ra = [Ra tc0*(trf.-tm).*itp_tca2d.(nearfloats.x, nearfloats.y)]
    ω = 2π/SOT.meanyear
    Ra = [Ra σannual.^2/2*cos.(ω*trf)/3.2 σannual.^2/2*sin.(ω*trf)/3.2] 
    Ra = [Ra σsemiannual.^2/2*cos.(2ω*trf)/3.2 σsemiannual.^2/2*sin.(2ω*trf)/3.2]

    Ea = [I(ma) cos.(ω*trf) sin.(ω*trf) cos.(2ω*trf) sin.(2ω*trf)]
    Ca = Ea*Ca*Ea'+N

    iCa = inv(Matrix(Ca))
    Pa = C0-Ra'*iCa*Ra
    xa = Ra'*iCa*nearfloats.Δsp1

    Rt = σs^2*exp.(-abs.(trd.-trde')/λt).*itp_c2d.(θd) 

    Rt = [Rt tc0*(trd.-tm).*itp_tc2d.(θd) zeros(m1,4)]
    Rt = [Rt; zeros(m1,m+5); zeros(4,m+1) R[end-3:end,end-3:end]]
    Rt = E1*Rt

    iCt = inv(Matrix(Ct))
    Pt = C0-Rt'*iCt*Rt
    xt = Rt'*iCt*[tpairs.Δτ; ppairs.Δτ]
    
    Rat = σs^2*exp.(-abs.(trf.-trd')/λt).*itp_ca.(nearfloats.x, nearfloats.y, θd')
    Rat += tc0*((trf.-tm)*(trd'.-tm)).*itp_tca.(nearfloats.x, nearfloats.y, θd')
    Rat = [Rat zeros(ma,m1) Ra[:,end-3:end]]
    Rat = Rat*E1'

    RiPat = (Ra'-Rt'*iCt*Rat')*inv(Ca-Rat*iCt*Rat')
    RiPta = (Rt'-Ra'*iCa*Rat)*inv(Ct-Rat'*iCa*Rat)
    xat = RiPat*nearfloats.Δsp1 + RiPta*[tpairs.Δτ; ppairs.Δτ]
    Pat = C0-(RiPat*Ra+RiPta*Rt)

    #iC = inv(Matrix([Ca Rat;Rat' Ct]))
    #Rat = [Ra;Rt]
    #Pat = C0-Rat'*iC*Rat
    #xat = Rat'*iC*[nearfloats.Δsp1; tpairs.Δτ; ppairs.Δτ]

    return Pa,xa,Pt,xt,Pat,xat,C0

end


"""
invert sub-basin average travel time anomalies using Argo and T waves
"""
function invertargotmt(grid,x,σtrend, σannual, σsemiannual,tm,nearfloats,tpairs,ppairs,tstation,pstations,evtpos,itp_c,itp_tc,itp_ca,itp_tca;tλx=200)

    σtrend /= meanyear
    tc0 = σtrend^2/itp_tc(0)
    λt,λx,σs,σm,σn,σp,σx,σ,σnp,σh = x

    N = σn^2*I

    # mean and trend
    m = size(grid,1)
    N = σn^2*I
    Ξ = Diagonal((σannual*ones(4)/3.2).^2/2)
    σm = σs
    
    xg,yg = grid[:,1],grid[:,2]
    dg = dist.(xg', yg', xg, yg)/1e3
    C0 = kron(Diagonal([σm^2, tc0]),exp.(-dg/tλx))

    trf = Dates.value.(nearfloats.t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
    ma = length(trf)
    dx = SOT.dist.(nearfloats.x', nearfloats.y', nearfloats.x, nearfloats.y)/1e3
    Ca = σs^2*sqrt(2)*exp.(-abs.(trf.-trf')/λt-dx/λx/sqrt(2))
    Ca = Ca.*cos.(dx/λx/sqrt(2).-π/4)
    Ca += (σm^2 .+ tc0*((trf.-tm)*(trf'.-tm))).*exp.(-dx/tλx)
    Ca = [Ca zeros(ma, 4); zeros(4, ma) Ξ]

    dga = SOT.dist.(xg', yg', nearfloats.x, nearfloats.y)/1e3
    Ra = σm^2*exp.(-dga/tλx)
    Ra = [Ra tc0*((trf.-tm)*ones(m)').*exp.(-dga/tλx)]
    ω = 2π/SOT.meanyear

    Ea = [I(ma) cos.(ω*trf) sin.(ω*trf) cos.(2ω*trf) sin.(2ω*trf)]
    Ca = Ea*Ca*Ea'+N

    iCa = inv(Matrix(Ca))
    Pa = C0-Ra'*iCa*Ra
    xa = Ra'*iCa*nearfloats.Δsp1
    
    ### twave
    _,_,θd = getazimuth(tpairs,ppairs,tstation[1],tstation[2],evtpos[1],evtpos[2])
    t1,E1,Ecs,Esn,Ecsh,Esnh = SOT.getE(tpairs,ppairs,tstation,pstations;hydro=true)

    nt,np = size(tpairs,1),size(ppairs,1)
    
    l = 1
    m1 = length(t1)
    # real time (days)
    trd = Dates.value.(t1 - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
    C = σs^2*exp.(-abs.(trd.-trd')/λt).*itp_c.(abs.((θd.-θd')))
    C += tc0*((trd.-tm)*(trd'.-tm)).*itp_tc.(abs.((θd.-θd')))
    E1 = [E1[:,1:2m1] E1[:,2m1+2:end]]

    # covariance matrix assuming no correlation between singular vector expansion coefficients
    R = [C zeros(l*m1, m1); zeros(m1, (l+1)*m1)] + σp^2*kron(sparse(ones(l+1, l+1)), I(m1))
    #R = [R zeros(size(R, 1), l); zeros(l, size(R, 2)) (σtrend.^2)]
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) I(2)*(σannual.^2/2)]
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) I(2)*(σsemiannual.^2/2)]

    # noise covariance
    Rcsn = Ecs*Ecs'+Esn*Esn'
    Rcsnh = Ecsh*Ecsh'+Esnh*Esnh'
    Nt = σ^2*spdiagm([ones(nt);zeros(np)])+σx^2*Rcsn+σh^2*Rcsnh
    Nt += σnp^2*spdiagm([zeros(nt);ones(np)])

    Ct = E1*R*E1'+Nt
    Rt = [zeros(m1,m) tc0*((trd.-tm)*ones(m)').*itp_tca.(xg', yg', θd)]
    Rt = E1*[Rt; zeros(m1+4,2m)]
    
    iCt = inv(Matrix(Ct))
    Pt = C0-Rt'*iCt*Rt
    xt = Rt'*iCt*[tpairs.Δτ; ppairs.Δτ]
    
    Rat = σs^2*exp.(-abs.(trf.-trd')/λt).*itp_ca.(nearfloats.x, nearfloats.y, θd')
    Rat += tc0*((trf.-tm)*(trd'.-tm)).*itp_tca.(nearfloats.x, nearfloats.y, θd')
    Rat = [Rat zeros(ma,m1) σannual.^2/2*cos.(ω*trf)/3.2 σannual.^2/2*sin.(ω*trf)/3.2] 
    Rat = [Rat σsemiannual.^2/2*cos.(2ω*trf)/3.2 σsemiannual.^2/2*sin.(2ω*trf)/3.2]
    Rat = Rat*E1'

    RiPat = (Ra'-Rt'*iCt*Rat')*inv(Ca-Rat*iCt*Rat')
    RiPta = (Rt'-Ra'*iCa*Rat)*inv(Ct-Rat'*iCa*Rat)
    xat = RiPat*nearfloats.Δsp1 + RiPta*[tpairs.Δτ; ppairs.Δτ]
    Pat = C0-(RiPat*Ra+RiPta*Rt)

    return Pa,xa,Pt,xt,Pat,xat

end
"""
invert 2d slowness anomalies using point data and T waves
"""
function invertargotag(grid,σtrend, σannual, tm,nearfloats,λt,λx,σs,σm,itp_tc;tλx=200)

    σtrend /= meanyear
    tc0 = σtrend^2/itp_tc(0)

    N = 0^2*I

    # mean and trend
    m = size(grid,1)

    Ξ = Diagonal((σannual*ones(4)/3.2).^2/2)
    #σm = σs
    
    xg,yg,tg = grid[:,1],grid[:,2],grid[:,3]
    dg = dist.(xg', yg', xg, yg)/1e3
    C0 = σs^2*sqrt(2)*exp.(-dg/λx/sqrt(2)).*cos.(dg/λx/sqrt(2).-π/4)
    #C0 += tc0*((tg.-tm)*(tg'.-tm)).*exp.(-dg/tλx)

    trf = Dates.value.(nearfloats.t - DateTime(2000, 1, 1))/1000/3600/24
    ma = length(trf)
    dx = SOT.dist.(nearfloats.x', nearfloats.y', nearfloats.x, nearfloats.y)/1e3
    Ca = σs^2*sqrt(2)*exp.(-abs.(trf.-trf')/λt-dx/λx/sqrt(2))
    Ca = Ca.*cos.(dx/λx/sqrt(2).-π/4)
    Ca += (σm^2 .+ tc0*((trf.-tm)*(trf'.-tm))).*exp.(-dx/tλx)
    Ca = [Ca zeros(ma, 4); zeros(4, ma) Ξ]

    dga = SOT.dist.(xg', yg', nearfloats.x, nearfloats.y)/1e3
    Ra = σs^2*sqrt(2)*exp.(-abs.(trf.-tg')/λt-dga/λx/sqrt(2)).*cos.(dga/λx/sqrt(2).-π/4)
    #Ra += tc0*((trf.-tm)*(tg'.-tm)).*exp.(-dga/tλx)
    ω = 2π/SOT.meanyear

    Ea = [I(ma) cos.(ω*trf) sin.(ω*trf) cos.(2ω*trf) sin.(2ω*trf)]
    Ca = Ea*Ca*Ea'+N

    iCa = inv(Matrix(Ca))
    Pa = C0-Ra'*iCa*Ra
    xa = Ra'*iCa*nearfloats.Δsp1

    return Pa,xa

end
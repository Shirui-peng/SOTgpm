include("../src/SOT.jl")
using .SOT, Printf, LinearAlgebra, Dates, HDF5, DataFrames, CSV, Statistics
using Random, Distributions, SparseArrays, PyPlot, Interpolations,NCDatasets
PyPlot.matplotlib[:rc]("mathtext",fontset="cm")        #computer modern font 
PyPlot.matplotlib[:rc]("font",family="STIXGeneral")
rc("font", size=8)
rc("axes", titlesize="medium")

# read argo data file
argofile = "results/argo/japan_1997-01_2022-12.csv"
floats = CSV.read(argofile, DataFrame)

# bulk sensitivity (s/K/km)
K0 = -6/3.2

# read rg grid
h5fl = "results/argo/kuroshio_rg.h5"
xrg = h5read(h5fl, "x")
yrg = h5read(h5fl, "y")
nxrg,nyrg=length(xrg),length(yrg)
nodes = (xrg, yrg)
h5fl2 = "results/argo/kuroshio_rg1922.h5"
trg = Millisecond.([h5read(h5fl, "t");h5read(h5fl2, "t")]) .+ DateTime(2000, 1, 1) 
Δsrg = cat(h5read(h5fl, "Δsp1"),h5read(h5fl2, "Δsp1");dims=3)
aidx = .!iszero.(Δsrg[:,:,1])
ntrg = length(trg)

# path travel time anomaly trend prior (s/yr)
σtrend = 0.037

# path travel time anomaly annual cycle prior (s)
σannual = 0.1

# path travel time anomaly semi-annual cycle prior (s)
σsemiannual = 0.1

# time scale (d)
λt = [60,60,57]

# length scale (km)
λx = [100/sqrt(2),100/sqrt(2),59]

# stochastic slowness anomaly scale (ms/km)
σs = [0.44,0.44,0.32]

# offset slowness anomaly scale (ms/km)
σm = [0.44,0.44]

# reference source location
evtpos = [142.85,38.10]

# source locations on path 
h11a0 = [[141.96,38.58],[142.06,38.87]]

# receiver locations
tstations = [[166.90986,19.71786],[166.652,19.283]]

# seismic data file name tails
tails = ["H11N3_x2y_16ex","WAKE_x2y_40ex"]

# receiver names
receivers = ["H11","WAKE"]

# grid number on path
nxk = 321

# dates for ssh data
sshs = ["20181124","20191123"]#,"20200318"]

# path-path trend prior covariance
θgrid = 0:0.1:20
cgrid = SOT.tcint.(θgrid)
itp_tc = linear_interpolation(θgrid, cgrid, extrapolation_bc=0)

# plot mapped temperature anomalies for rg, this study, and ssh
fig,ax=subplots(2,3,figsize=(190/25.4, 190/25.4-1.7),sharex=true,sharey=true)
#fig.suptitle("Feb. 15, 2005 (upper row) vs. Nov 24, 2018 (lower row)")#"vs. Mar. 18, 2020 (lower row)")
fig.text(0.46,0.9,"Nov. 24, 2018")
fig.text(0.46,0.505,"Nov. 23, 2019")
for (i,te) in enumerate([DateTime(2018, 11, 24),DateTime(2019, 11, 23)])#,DateTime(2020, 3, 18)])
    
    # labels
    ax[i,1].set_title("Temperature anomaly (RG)")
    ax[i,2].set_title("Temperature anomaly (this study)")
    ax[i,3].set_title("Sea level anomaly")
    
    # events to invert
    m = 1
    lone,late,azme = ones(m)*h11a0[i][1],ones(m)*h11a0[i][2],zeros(m)
    nms = ["longitude","latitude","time","azimuth"]
    events = DataFrame([lone late [te] azme], nms)

    # select nearby floats
    Δt = Day(2*λt[i])
    t0,t1,x0,y0,z0 = te-Δt,te+Δt,141,36,-1.9e3
    #floats = floats[(floats.z0.<z0),:]
    ifloat = ((floats.y .<=y0) .| (floats.x .>= x0)).&(t0 .≤ floats.t .< t1).& (abs.(floats.Δsp1) .< 10)
    nearfloats = unique(floats[ifloat,:], [:x,:y,:t])

    # location of max anomaly
    nt = findfirst(trg .>= te)
    ijm = Tuple(argmax(-Δsrg[xrg.<150,:,nt]))
    println(ijm)
    println([xrg[ijm[1]],yrg[ijm[2]]])

    # read ssh
    filename = @sprintf("data/ssh/dt_global_twosat_phy_l4_%s_vDT2021.nc",sshs[i])
    fl = Dataset(filename, "r")
    lon = Array{Float64}(fl["longitude"][:])
    lat = Array{Float64}(fl["latitude"][:])
    ix,iy = 135 .<= lon .<= 170,15 .<= lat .<= 45
    sla = Array{Float64}(replace(fl["sla"][ix,iy,1],(missing=>NaN)))
    mask = replace(isnan.(sla),(false=>NaN))

    # path grid
    xk0, yk0 = SOT.findpath(tstations[i], h11a0[i], nxk)
    
    # rg map
    im1=ax[i,1].pcolormesh(xrg,yrg,Δsrg[:,:,nt]'/K0,vmin=-1,vmax=1,
                         cmap="RdBu_r",shading="auto",rasterized=true)
    
    # plot boundaries and float locations
    for j = 1:3
      ax[i,j].plot(xk0,yk0,lw=1,c="k",zorder=3)
      ax[i,j].pcolormesh(lon[ix],lat[iy],mask',cmap="Greys",vmin=0,vmax=2,
                       shading="auto",rasterized=true)
      ax[i,j].set_xlim([135,170])
      Rt = exp.(-abs.(Dates.value.(nearfloats.t-te)/1e3/24/3600)/λt[i])
      ax[i,j].scatter(nearfloats.x,nearfloats.y, s=3,c="k",
                      alpha=Rt/max(Rt...),zorder=3, linewidths=0)
      ax[i,j].set_title("($(('a':'z')[3*(i-1)+j]))", loc="left")
      ax[i,j].set_aspect(1/cosd(30))
      ax[end,j].set_xlabel("longitude")
    end
    
    # grid to invert for
    ted = Dates.value.(te - DateTime(2000, 1, 1))/1000/3600/24
    gxy = zeros(sum(aidx),3)
    for j = 1:nxrg, k = 1:nyrg
        if aidx[j,k] 
            na = 1 + sum(gxy[:,3].>0)
            gxy[na,:]=[xrg[j],yrg[k],ted]
        end
    end
    
    # invert for stochastic temperature anomaly map
    ttw = h5read(@sprintf("results/anomalies/japan_%s.h5",tails[i]), "t")
    tm = (min(ttw...)+max(ttw...))/2/1000/3600/24
    Pa,xa = SOT.invertargotag(gxy,3.7σtrend,σannual,tm,nearfloats,λt[i],λx[i],σs[i],σm[i],itp_tc)
    
    # assign inverted results
    gae = fill(NaN, nxrg, nyrg, 3,2)
    ng=0
    for i = 1:nxrg, j = 1:nyrg
        if aidx[i,j]
            ng+=1
            gae[i,j,1,1],gae[i,j,1,2]=xa[ng],sqrt(diag(Pa)[ng])
        end
    end

    # add trend
    h5name = @sprintf("results/argo/kuroshio_mt_%s.h5",receivers[i])
    gtd = h5read(h5name, "gtd")
    DT = (ted-tm)*gtd[:,:,1]'/K0
    
    # plot temperature anomaly and ssh maps
    im2=ax[i,2].pcolormesh(xrg,yrg,gae[:,:,1,1]'/K0 .+ DT,vmin=-1,vmax=1,
                          cmap="RdBu_r",shading="auto",rasterized=true)
    im3=ax[i,3].pcolormesh(lon[ix],lat[iy],sla',vmin=-1,vmax=1,
                          cmap="RdBu_r",shading="auto",rasterized=true)

    # add texts, colorbars, and sector areas
    ax[i,1].set_ylabel("latitude")
    if i==1
        ax[1,1].text(xk0[1]-2.5,yk0[1]-2,"Wake")
        ax[1,1].text(xk0[1]-2.5,yk0[1]-4,"Island")
        ax[1,1].text(xk0[end]-6.5,yk0[end]+2,"Japan")
        y1,dy = .58,-1.5
        cbaxes = fig.add_axes([0.14, y1, 0.06, 0.02]) 
        cbaxes.text(1.3, dy,"(K)")
        cbar = fig.colorbar(im1, cax=cbaxes, ticks=[-1,0,1], orientation="horizontal")
        cbaxes = fig.add_axes([0.40, y1, 0.06, 0.02]) 
        cbaxes.text(1.3, dy,"(K)")
        cbar = fig.colorbar(im2, cax=cbaxes, ticks=[-1,0,1], orientation="horizontal")
        cbaxes = fig.add_axes([0.66, y1, 0.06, 0.02]) 
        cbaxes.text(1.3, dy,"(m)")
        cbar = fig.colorbar(im3, cax=cbaxes, ticks=[-1,0,1], orientation="horizontal")
    end
    a0 = SOT.azimuth.(tstations[i][2],tstations[i][1],evtpos[2],evtpos[1])
    arg = SOT.azimuth.(tstations[i][2],tstations[i][1],yrg',xrg) .- a0
    drg = SOT.dist.(tstations[i][1],tstations[i][2], xrg, yrg')/1e6
    ai = -6.5 .<= arg .<= 3.5
    ai = ai .& (drg.<=3.2) .& (cosd.(arg).>=0)
    gxy = zeros(sum(ai),2)
    for n = 1:nxrg, m = 1:nyrg
        if ai[n,m] 
            na = 1 + sum(gxy[:,2].>0)
            gxy[na,:]=[xrg[n],yrg[m]]
        end
    end
    iW,iN,iE,iS = argmin(gxy[:,1]),argmax(gxy[:,2]),argmax(gxy[:,1]),argmin(gxy[:,2])
    yp = [gxy[iW,2],gxy[iN,2],gxy[iE,2],gxy[iS,2],gxy[iW,2]]
    xp = [gxy[iW,1],gxy[iN,1],gxy[iE,1],gxy[iS,1],gxy[iW,1]]
    xk1,yk1 = SOT.findpath([xp[2],yp[2]],tstations[i],nxk)
    xk2,yk2 = SOT.findpath([xp[1],yp[1]],tstations[i],nxk)
    ax[i,1].plot(xk1,yk1,c="k",lw=.8)
    ax[i,1].plot(xk2,yk2,c="k",lw=.8)
    ax[i,1].plot(xp[1:2],yp[1:2],c="k",lw=.8)

end
subplots_adjust(top=0.9,hspace=0.0,wspace=0.05)
fig.savefig(@sprintf("results/argo/kuroshio_rg_argo_ssh%d.pdf",length(sshs)),bbox_inches="tight",dpi=300)
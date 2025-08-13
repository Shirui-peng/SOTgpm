include("../src/SOT.jl")
using .SOT, PyPlot, Printf, Dates, LinearAlgebra, Statistics, SparseArrays, DataFrames, CSV
using HDF5, Interpolations, Random, Distributions, NCDatasets, StatsBase

# identifier for experiment
eqname = "japan"
tsname = "H11ccx"
evtpos = [38.10, 142.85]

# P-wave (reference) stations
pstations = ["IU.MAJO.00.BHZ","IU.MAJO.10.BHZ","PS.TSK..BHZ", "II.ERM.00.BHZ", "G.INU.00.BHZ"] 
pstnlats,pstnlons = [36.55,36.55,36.21,42.02,35.35],[138.2,138.2,140.11,143.16,137.03]
pstalocs = DataFrame(station=pstations,slat=pstnlats,slon=pstnlons)

# P-wave download source
psrc = "IRIS"

# intervals to which to cut P waveforms
pintervals = [[-3, 47], [-3, 47], [-3, 47], [-3, 47], [-3, 47]]

# frequency bands to which to filter P waveforms
pfreqbands = [[1.5, 2.5], [1.5, 2.5], [1, 3], [1.5, 2.5], [1.5, 2.5]]

# T-wave station
tstations = ["H11N3"]
tstalocs = [19.71786, 166.90986] 

# T-wave download source
tsrc = "IRIS"

# T-wave time window around predicted arrival time
tintervals = [[-20, 40]]

# T-wave filtering window width
tavgwidth = 0.5

# T-wave reference frequency at which to find first max CC
treffreq = 2.5

# frequencies used in inversion
tinvfreq = [2.5, 3.5]

# minimum CCs for T-wave pairs (at inversion frequencies)
tmincc = [0.6, 0.3]

# download P-wave data
#SOT.downloadseisdata(eqname, tsname, pstations; src=psrc)

# cut and filter P waveforms
#SOT.cutpwaves(eqname, tsname, pstations, pintervals, pfreqbands)

# find P-wave pairs
SOT.findpairs(eqname, tsname, pstations, pintervals, pfreqbands; saveplot=true)#; kstart=51748)

# download T-wave data
#SOT.downloadseisdata(eqname, tsname, ["IM.H11N3..EDH"]; src=tsrc)
 
# measure T-wave lags Δτ
#SOT.twavepickold(eqname, tsname, tstations, tintervals, tavgwidth, treffreq, pstations, pintervals, pfreqbands;soundspeed=1.47e3)

# manually exclude pairs
#excludepairs = CSV.read("data/catalogs/japan_H11_exclude.csv", DataFrame)

# collect usable pairs
#tpairs, ppairs = SOT.collectpairs(eqname, tsname, tstations, tintervals, tavgwidth, treffreq,
#                                  tinvfreq, tmincc, pstations, pintervals, pfreqbands;soundspeed=1.47e3)#,excludepairs=excludepairs)

ppairs.stn6 = [s[1:6] for s in ppairs.station]
ppairs = unique(ppairs, [:stn6,:event1,:event2])
#select!(ppairs, Not(:stn6))

CSV.write(@sprintf("results/pairs/japan_%s_tpairs_%.1fa%.1fhz_full.csv",tsname,tinvfreq[1],tinvfreq[2]), tpairs)
CSV.write(@sprintf("results/pairs/japan_%s_ppairs_%.1fa%.1fhz_full.csv",tsname,tinvfreq[1],tinvfreq[2]), ppairs)

excludeevts = false
if excludeevts
    # load event catalog
    events = DataFrame(CSV.File(@sprintf("data/catalogs/%s_%s.csv", eqname, tsname)))
    # load and combine catalogs of P-wave pairs
    eevents = DataFrame[]
    for e = eachrow(excludepairs)
      exclude = (events.time .== e.event1) 
      push!(eevents, events[exclude,:])
      exclude = (events.time .== e.event2) 
      push!(eevents, events[exclude,:])
    end
    eevents = sort(unique(vcat(eevents...)))
    CSV.write("data/catalogs/japan36s_H11.csv", eevents)
end

# number of frequencies
l = length(tinvfreq)-1

# correlation time (days)
#λt = 60

# correlation azimuth (degrees)
λθ = 2.0

# solution standard deviation for travel time anomalies (s)
στ = [0.27]

# location noise (s)
σx,σh = 0.019,0.028

# noise (s)
σn,σnp = 8.1e-3,2.2e-3

# origin time correction standard deviation (s)
σp = 0.974

# trend prior for coefficients of singular vectors (s/day)
σtrend = 0.01/SOT.meanyear

# annual cycle prior (s)
σannual = 0.1

# semi-annual cycle prior (s)
σsemiannual = 0.1

# trend prior for coefficients of singular vectors (s/day)
#σθtrend = 0.01
@printf("lt = %.0f days, n = %.2e s, nx = %.2e s\n\n",λt,σn,σx)

t, lon, lat, θ, E, R, N, P, D, invRyy = SOT.invertf1(tpairs, ppairs, tstalocs, pstalocs, evtpos, λt, λθ, στ, σx, σh, σn, σnp, σp; σtrend, σannual, σsemiannual)

# number of good T- and P-wave pairs
nt = size(tpairs, 1)
np = size(ppairs, 1)
m = length(t)
    
tpairs.Δτ,tpairs.Δτp,tpairs.cs,tpairs.x1,tpairs.x2,tpairs.Δτ2,zint,zinp,zio = SOT.correctcycleskippingf1r2(eqname, tstations, tpairs, ppairs, E, R, N, P, m)

qplot = true
if qplot
    p = 0.5:0.5:99.5
    α = 0.05
    qs1,qs2,qs3 = percentile(zio,p),percentile(zint,p),percentile(zinp,p)
    qt = quantile.(Normal(), p/100)
    lb,ub = zeros(length(p),3),zeros(length(p),3)
    for i = 1:length(p)
        for (j,ns) in enumerate([m, nt, np])
            d = NoncentralT(ns-1,-sqrt(ns)*qt[i])
            lb[i,j],ub[i,j] = -quantile(d, 1-α/2)/sqrt(ns),-quantile(d, α/2)/sqrt(ns)
        end
    end
    # small font
    rc("font", size=8)
    rc("axes", titlesize="medium")
    fig,ax=plt.subplots(figsize=(4.8,3.6))
    ax.plot(qt,qs1-qt,label="origin time error",c="tab:blue")
    ax.fill_between(qt, (lb[:,1]-qt), (ub[:,1]-qt), alpha=.2, zorder=3, color="tab:blue", linewidth=0)
    ax.plot(qt,qs2-qt,label="\$T\$-wave residual",c="tab:orange")
    ax.fill_between(qt, (lb[:,2]-qt), (ub[:,2]-qt), alpha=.2, zorder=3, color="tab:orange", linewidth=0)
    ax.plot(qt,qs3-qt,label="\$P\$-wave residual",c="tab:green")
    ax.fill_between(qt, (lb[:,3]-qt), (ub[:,3]-qt), alpha=.2, zorder=3, color="tab:green", linewidth=0)
    ax.axhline(0,color="black",ls=":",lw=1)
    ax.set_xlabel("theoretical quantile")
    ax.set_ylabel("sample quantile - theoretical quantile")
    ax.legend(frameon=false)
    fig.tight_layout()
    fig.savefig(@sprintf("results/plots/quantile_%s_%s.pdf", eqname, tsname))
end

# collect delays into data vector
y = [reshape(vcat([(tpairs.Δτ[i])' for i = 1:nt]...), l*nt); ppairs.Δτ]
y2 = [tpairs.Δτ2; ppairs.Δτ]

tr = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24

# invert
a = P*E'*inv(Array(N))*y
a2 = P*E'*inv(Array(N))*y2
 
# extract trends
trends = a[2*m+1]
Ptrends = P[2*m+1,2*m+1]
@printf("\ntrend = %.2e K/yr, ptrend = %.2e\n",-trends*SOT.meanyear/6,2*sqrt(Ptrends)*SOT.meanyear/6)

# extract trends
annual = a[2*m+2:2*m+3]
Pannual = diag(P[2*m+2:2*m+3,2*m+2:2*m+3])
@printf("\nannual = %s K, pannual = %s K\n",-annual/6,2*sqrt.(Pannual)/6)

# extract trends
semiannual = a[2*m+4:2*m+5]
Psannual = diag(P[2*m+4:end,2*m+4:end])
@printf("\nsemiannual = %s K, psannual = %s K\n\n",-semiannual/6,2*sqrt.(Psannual)/6)

# extract azimuth trends
#atrends = a[end]
#Patrends = P[end,end]
#@printf("\natrend = %.2e s/deg, patrend = %.2e\n",atrends,2*sqrt(Patrends))

#annual = [-0.03396737591858773, 0.053609609834979394] s, pannual = [0.05210037025550977, 0.052217578713336724] s
#semiannual = [-0.033771452928324903, -0.023933543776224326] s, psannual = [0.037785276755231814, 0.039253873902431166] s

# plot residuals with and without smoothing
n1 = y - E*a
n2 = y - E*(E\y)
tpairs.n1,tpairs.n2 = n1[1:nt],n2[1:nt]
ppairs.n1,ppairs.n2 = n1[nt+1:nt+np],n2[nt+1:nt+np]
n12,n22 = y2 - E*a2,y2 - E*(E\y2)

tn1maxi = argmax(abs.(n12[1:nt]))
@printf("tpair n12 std: %.3f; outlier: %s %s, Dt=%s, Dtp=%.3fs, n1=%.2f s\n", std(n12), tpairs.event1[tn1maxi],tpairs.event2[tn1maxi],tpairs.Δτ2[tn1maxi],tpairs.Δτp[tn1maxi],n12[tn1maxi])
tn2maxi = argmax(abs.(n22[1:nt]))
@printf("tpair n22 std: %.3f; outlier: %s %s, Dt=%s, Dtp=%.3fs, n2=%.2f s\n\n", std(n22), tpairs.event1[tn2maxi],tpairs.event2[tn2maxi],tpairs.Δτ2[tn2maxi],tpairs.Δτp[tn2maxi],n22[tn2maxi])

tn1maxi = argmax(abs.(tpairs.n1))
@printf("tpair n1 std: %.3f; outlier: %s %s, Dt=%s, Dtp=%.3fs, cs=%d, n1=%.2f s\n", std(tpairs.n1), tpairs.event1[tn1maxi],tpairs.event2[tn1maxi],tpairs.Δτ[tn1maxi],tpairs.Δτp[tn1maxi],tpairs.cs[tn1maxi],tpairs.n1[tn1maxi])
tn2maxi = argmax(abs.(tpairs.n2))
@printf("tpair n2 std: %.3f; outlier: %s %s, Dt=%s, Dtp=%.3fs, n2=%.2f s\n\n", std(tpairs.n2), tpairs.event1[tn2maxi],tpairs.event2[tn2maxi],tpairs.Δτ[tn2maxi],tpairs.Δτp[tn2maxi],tpairs.n2[tn2maxi])
pn1maxi = argmax(abs.(ppairs.n1))
@printf("ppair n1 std: %.3f; outlier: %s %s, n1=%.2f s\n", std(ppairs.n1), ppairs.event1[pn1maxi],ppairs.event2[pn1maxi],ppairs.n1[pn1maxi])
pn2maxi = argmax(abs.(ppairs.n2))
@printf("ppair n2 std: %.3f; outlier: %s %s, n2=%.2f s\n", std(ppairs.n2), ppairs.event1[pn2maxi],ppairs.event2[pn2maxi],ppairs.n2[pn2maxi])
if abs(ppairs.n2[pn2maxi])<0.1
  if abs(tpairs.n1[tn1maxi]) > abs(tpairs.n2[tn2maxi])
    ppairsexclude = ppairs[(ppairs.event1 .== tpairs.event1[tn1maxi]) .&& (ppairs.event2 .== tpairs.event2[tn1maxi]), :]
  else
    ppairsexclude = ppairs[(ppairs.event1 .== tpairs.event1[tn2maxi]) .&& (ppairs.event2 .== tpairs.event2[tn2maxi]), :]
  end
else
  ppairsexclude = ppairs[(ppairs.event1 .== ppairs.event1[pn2maxi]) .&& (ppairs.event2 .== ppairs.event2[pn2maxi]), :]
end
ppairsexclude = ppairsexclude[:, [:station, :event1, :event2]]
tpairs.Δτ = vec(vcat([(tpairs.Δτ[i])' for i = 1:nt]...)) 
CSV.write(@sprintf("results/pairs/japan_%s_ppairsexclude.csv",tstations[1]), ppairsexclude)
CSV.write(@sprintf("results/pairs/japan_%s_ppairs_2.5a3.5hz.csv",tstations[1]), ppairs)
CSV.write(@sprintf("results/pairs/japan_%s_tpairs_2.5a3.5hz.csv",tstations[1]), tpairs)

hist = true
if hist
    b = -0.2:5e-3:0.2
    @assert all((b[1] .< n1 .< b[end]) .& (b[1] .< n2 .< b[end]))
    fig, ax = subplots(2, l+2; sharex=true, sharey=true)
    for i = 1:l
      ax[1,i].hist(n1[(i-1)*nt+1:i*nt]; bins=b)
      ax[2,i].hist(n2[(i-1)*nt+1:i*nt]; bins=b)
      ax[1,i].set_title("T wave, freq. $i")
      ax[2,i].set_title("T wave, freq. $i")
    end
    ax[1,l+1].hist(n1[l*nt+1:l*nt+np]; bins=b)
    ax[2,l+1].hist(n2[l*nt+1:l*nt+np]; bins=b)
    ax[1,l+1].set_title("P wave")
    ax[2,l+1].set_title("P wave")
    ax[1,1].set_yscale("log")
    ax[1,1].set_ylabel("count")
    ax[2,1].set_ylabel("count")
    
    ax[1,3].hist(n12[1:nt]; bins=b)
    ax[2,3].hist(n22[1:nt]; bins=b)
    ax[1,3].set_title("T wave, f2")
    ax[2,3].set_title("T wave, f2")
    for i = 1:l+1
      ax[2,i].set_xlabel("residual (s)")
    end
    fig.tight_layout()
    fig.savefig(@sprintf("results/plots/japan_%s_residuals_hist.pdf",tstations[1]))
    @printf("histogram finished\n")
end

# reconstruct full travel time anomalies
τ = reshape(D*a, (m, 1)) 
τ2 = reshape(D*(E\y), (m, 1))
eτ = reshape(sqrt.(diag(D*P*D')), (m, 1))

filename = @sprintf("results/pairs/japan_%s_x2y.h5",tstations[1])
# save to file
h5open(filename, "w") do file
  write(file, "t", Dates.value.(t .- DateTime(2000, 1, 1, 0, 0, 0)))
  write(file, "lon", lon)
  write(file, "lat", lat)
  write(file, "θ", θ)
  write(file, "τ", τ)
  write(file, "e", eτ)
  write(file, "y", y)
  write(file, "a", a)
end

df_T = DataFrame(time = vec(t), azimuth = vec(θ), τ = vec(τ), τ2 = vec(τ2), e = vec(eτ))
CSV.write(@sprintf("results/anomalies/japan_%s.csv",tstations[1]), df_T)

#td = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24 

# interpolate onto regular grid
θi = -9:0.25:9 
θd = θ

# interpolate onto regular grid
plot = false
for year = 2008:2021
  @printf("\ninterpolation year %d\n",year)
  ti = DateTime(year, 1, 1, 12, 0, 0) : Day(10) : DateTime(year, 12, 31, 12, 0, 0)
  τi,ei,ei0 = SOT.regulargrid(td, θd, ti, θi, a, R, invRyy, E, λt, λθ, στ)
  if year == 2008
    global tiall,τiall,eiall,ei0all = ti,τi,ei,ei0
  else    
    tiall,τiall,eiall,ei0all = vcat(tiall,ti),hcat(τiall,τi),hcat(eiall,ei),hcat(ei0all,ei0)
    @printf("ti shape: %s\n",size(tiall))
  end
  if plot
    ri = eiall./στ
    ri[ri .> 0.5] .= 1
    ri[ri .<= 0.5] .= NaN
    fig, ax = subplots(2,1,figsize=(190/25.4, 190/25.4),sharex=true)
    ax[1].pcolormesh(tiall,θi,τiall,vmin=-1,vmax=1,cmap="RdBu",shading="nearest",zorder=0)
    ax[1].pcolormesh(tiall,θi,ri,vmin=0,vmax=1.5,cmap="Greys",alpha=0.2,shading="nearest")
    ax[2].pcolormesh(tiall,θi,eiall./στ,vmin=0,vmax=1,shading="nearest")
    for i = 1:2
      ax[i].scatter(df_T.time,df_T.azimuth,s=5,color="black")
      ax[i].set_ylabel("azimuth (deg)")
      ax[i].set_xlim([tiall[1],tiall[end]])
    end
    fig.tight_layout()
    fig.savefig(@sprintf("results/plots/japan_%s_regularmap_%d.pdf",tstations[1],year))
  end
end

tiall = Dates.value.(tiall - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24

filename = @sprintf("results/anomalies/japan_%s_map.nc",tstations[1])
ds = Dataset(filename,"c")

# Define the dimension 
defDim(ds,"t",length(tiall))
defDim(ds,"a",length(θi))

# Define the variables 
vt = defVar(ds,"time",Float64,("t",))
va = defVar(ds,"azimuth",Float64,("a",)) 
vτ = defVar(ds,"tau",Float64,("a","t"))
ve = defVar(ds,"e",Float64,("a","t"))
vre = defVar(ds,"re",Float64,("a","t"))
vt[:] = tiall
va[:] = θi
vτ[:,:] = τiall
ve[:,:] = eiall
vre[:,:] = eiall./ei0all
close(ds)

sp = false
if sp
    d = Normal()
    Random.seed!(123) # Setting the seed
    
    r(n) = 2/π*asin.(sqrt.(((1:n) .- 0.5)./n))
    s(d,y,μ,σ) = 2/π*asin.(sqrt.(cdf.(d,(y .- μ)./σ)))
    
    fig, ax = subplots(2, l+1; sharex=true, sharey=true)
    for i = 1:l+1, j = 1:l+1
      ax[i,j].plot([0,1],[0,1],color="k")
    end
    xi = r(nt)
    yref = s(d,sort(rand(d, nt)),0,1)
    ax[1,1].scatter(xi,yref,s=3,color="r")
    ax[2,1].scatter(xi,yref,s=3,color="r")
    sort!(tpairs,[:n1])
    yi = s(d,tpairs.n1,0,σn)
    ax[1,1].scatter(xi,yi,s=3)
    sort!(tpairs,[:n2])
    yi = s(d,tpairs.n2,0,sqrt(sum(tpairs.n2.^2)/nt))
    ax[2,1].scatter(xi,yi,s=3)
    
    xi = r(np)
    yref = s(d,sort(rand(d, np)),0,1)
    ax[1,2].scatter(xi,yref,s=3,color="r")
    ax[2,2].scatter(xi,yref,s=3,color="r")
    sort!(ppairs,[:n1])
    yi = s(d,ppairs.n1,0,σn)
    ax[1,2].scatter(xi,yi,s=3)
    sort!(ppairs,[:n2])
    yi = s(d,ppairs.n2,0,sqrt(sum(ppairs.n2.^2)/np))
    ax[2,2].scatter(xi,yi,s=3)
    fig.tight_layout()
    fig.savefig(@sprintf("results/plots/japan_%s_residuals_sp.pdf",tstations[1]))
end

autocor = false
if autocor
    fig, ax = subplots(2,1)
    autocor(n,lags) = [sum(n[1:length(n) - abs(lags[i])].*n[1 + abs(lags[i]):end])/sum(n.^2) for i=1:length(lags)]
    lags = -(nt-10):(nt-10)
    ax[1].plot(lags,autocor(tpairs.n1, lags))
    ax[1].plot(lags,autocor(tpairs.n2, lags))
    lags = -(np-10):(np-10)
    ax[2].plot(lags,autocor(ppairs.n1, lags))
    ax[2].plot(lags,autocor(ppairs.n2, lags))
    fig.tight_layout()
    fig.savefig(@sprintf("results/plots/japan_%s_residuals_autocor.pdf",tstations[1]))
end
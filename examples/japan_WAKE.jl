include("../src/SOT.jl")
using .SOT, PyPlot, Printf, Dates, LinearAlgebra, Statistics, SparseArrays, DataFrames, CSV
using HDF5, Interpolations, Random, Distributions,NCDatasets

# identifier for experiment
eqname = "japan"  
tsname = "WAKE"
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
tstations = ["IU.WAKE.00.BHZ","IU.WAKE.10.BHZ","IU.WAKE.20.BHZ"]
tstalocs = [19.283, 166.652]

# T-wave download source
tsrc = "IRIS"

# T-wave time window around predicted arrival time
tintervals = [[-15, 25],[-15, 25],[-15, 25]]

# T-wave filtering window width
tavgwidth = 0.5

# T-wave reference frequency at which to find first max CC
treffreq = 2.5

# frequencies used in inversion
tinvfreq = [2.5, 3.5]

# minimum CCs for T-wave pairs (at inversion frequencies)
tmincc = [0.6, 0.3]

# maximum separations in lat/lon/mag for pairs
# Δllm = [2.0,2.5,1.5]

# download P-wave data
#SOT.downloadseisdata(eqname, tsname, pstations; src=psrc)

# cut and filter P waveforms
#SOT.cutpwaves(eqname, tsname, pstations, pintervals, pfreqbands)
 
# find P-wave pairs
#SOT.findpairs(eqname, tsname, pstations, pintervals, pfreqbands; kstart=27654)

# download T-wave data
#SOT.downloadseisdata(eqname, tsname, tstations; src=tsrc)

# measure T-wave lags Δτ
#SOT.twavepickold(eqname, tsname, tstations, tintervals, tavgwidth, treffreq, pstations, pintervals, pfreqbands; soundspeed=1.47e3,saveplot=true)

#pn2max = 3.0
#while abs(pn2max) > 0.05
# manually exclude pairs
#excludepairs = CSV.read("data/catalogs/japan_WAKE_exclude.csv", DataFrame)
# collect usable pairs
tpairs, ppairs = SOT.collectpairs(eqname, tsname, tstations, tintervals, tavgwidth, treffreq, tinvfreq, tmincc, pstations, pintervals, pfreqbands;soundspeed=1.47e3,l2=true)#,excludepairs=excludepairs)

tpairs = unique(tpairs, [:event1,:event2])
ppairs.stn6 = [s[1:6] for s in ppairs.station]
#ppairsu = unique(ppairs, [:stn6,:event1,:event2])
select!(ppairsu, Not(:stn6))

#CSV.write(@sprintf("results/pairs/japan_%s_tpairs_%.1fa%.1fhz_full.csv",tsname,tinvfreq[1],tinvfreq[2]), tpairs)

# number of frequencies
l = length(tinvfreq)-1

# correlation time (days)
λt = 59

# correlation azimuth (degrees)
λθ = 2.2

# solution standard deviation for travel time anomalies (s)
στ = [0.29]

# location noise (s)
σx,σh = 0.028,0

# noise (s)
σn,σnp = 10e-3,2.9e-3

# origin time correction standard deviation (s)
σp = 0.86

# trend priors for coefficients of singular vectors (s/day)
σtrend = 0.01/SOT.meanyear

# annual cycle prior (s)
σannual = 0.1

# semi-annual cycle prior (s)
σsemiannual = 0.1

# trend prior for coefficients of singular vectors (s/day)
#σθtrend = 0.01
@printf("lt = %.0f days, n = %.2e s, nx = %.2e s\n\n",λt,σn,σx)

t, lon, lat, θ, E, R, N, P, D, invRyy = SOT.invertf1(tpairs, ppairsu, tstalocs, pstalocs, evtpos, λt, λθ, στ, σx, σh, σn, σnp, σp; σtrend, σannual, σsemiannual)

# number of good T- and P-wave pairs
nt = size(tpairs, 1)
np = size(ppairsu, 1)
m = length(t) 
    
tpairs.Δτ,tpairs.Δτp,tpairs.cs,tpairs.x1,tpairs.x2,tpairs.Δτ2 = SOT.correctcycleskippingf1r2(eqname, tstations, tpairs, ppairsu, E, R, N, P, m; l2=true)

# collect delays into data vector
y = [reshape(vcat([(tpairs.Δτ[i])' for i = 1:nt]...), l*nt); ppairsu.Δτ]
y2 = [tpairs.Δτ2; ppairsu.Δτ]

tr = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24

# invert
a = P*E'*inv(Array(N))*y
a2 = P*E'*inv(Array(N))*y2
 
# extract trends
trends = a[2*m+1]
Ptrends = P[2*m+1,2*m+1]
@printf("\ntrend = %.2e K/yr, ptrend = %.2e K/yr\n",-trends*SOT.meanyear/6,2*sqrt(Ptrends)*SOT.meanyear/6)

# extract trends
annual = a[2*m+2:2*m+3]
Pannual = diag(P[2*m+2:2*m+3,2*m+2:2*m+3])
@printf("\nannual = %s K, pannual = %s K\n",annual/6,2*sqrt.(Pannual)/6)

# extract trends
semiannual = a[2*m+4:2*m+5]
Psannual = diag(P[2*m+4:end,2*m+4:end])
@printf("\nsemiannual = %s K, psannual = %s K\n\n",semiannual/6,2*sqrt.(Psannual)/6)

# extract azimuth trends
#atrends = a[end]
#Patrends = P[end,end]
#@printf("\natrend = %.2e s/deg, patrend = %.2e\n",atrends,2*sqrt(Patrends))

#annual = [-0.01614593182503743, 0.031375631554657835] s, pannual = [0.05219946809286837, 0.05301471504447643] s
#semiannual = [-0.037620439290754426, 0.00018506178054533538] s, psannual = [0.04126444113991889, 0.040768364899831565] s

# plot residuals with and without smoothing
n1 = y - E*a
n2 = y - E*(E\y)
tpairs.n1,tpairs.n2 = n1[1:nt],n2[1:nt]
ppairsu.n1,ppairsu.n2 = n1[nt+1:nt+np],n2[nt+1:nt+np]
n12,n22 = y2 - E*a2,y2 - E*(E\y2)

tn1maxi = argmax(abs.(n12[1:nt]))
@printf("tpair n12 std: %.3f; outlier: %s %s, Dt=%s, Dtp=%.3fs, n1=%.2f s\n", std(n12), tpairs.event1[tn1maxi],tpairs.event2[tn1maxi],tpairs.Δτ2[tn1maxi],tpairs.Δτp[tn1maxi],n12[tn1maxi])
tn2maxi = argmax(abs.(n22[1:nt]))
@printf("tpair n22 std: %.3f; outlier: %s %s, Dt=%s, Dtp=%.3fs, n2=%.2f s\n\n", std(n22), tpairs.event1[tn2maxi],tpairs.event2[tn2maxi],tpairs.Δτ2[tn2maxi],tpairs.Δτp[tn2maxi],n22[tn2maxi])

tn1maxi = argmax(abs.(tpairs.n1))
@printf("tpair n1 std: %.3f; outlier: %s %s, Dt=%s, Dtp=%.3fs, n1=%.2f s\n", std(tpairs.n1), tpairs.event1[tn1maxi],tpairs.event2[tn1maxi],tpairs.Δτ[tn1maxi],tpairs.Δτp[tn1maxi],tpairs.n1[tn1maxi])
tn2maxi = argmax(abs.(tpairs.n2))
@printf("tpair n2 std: %.3f; outlier: %s %s, Dt=%s, Dtp=%.3fs, n2=%.2f s\n\n", std(tpairs.n2), tpairs.event1[tn2maxi],tpairs.event2[tn2maxi],tpairs.Δτ[tn2maxi],tpairs.Δτp[tn2maxi],tpairs.n2[tn2maxi])
pn1maxi = argmax(abs.(ppairsu.n1))
@printf("ppair n1 std: %.3f; outlier: %s %s, n1=%.2f s\n", std(ppairsu.n1), ppairsu.event1[pn1maxi],ppairsu.event2[pn1maxi],ppairsu.n1[pn1maxi])
pn2maxi = argmax(abs.(ppairsu.n2))
@printf("ppair n2 std: %.3f; outlier: %s %s, n2=%.2f s\n", std(ppairsu.n2), ppairsu.event1[pn2maxi],ppairsu.event2[pn2maxi],ppairsu.n2[pn2maxi])
if abs(ppairsu.n2[pn2maxi])<0.1
  if abs(tpairs.n2[tn2maxi]) < 0.1
    ppairsexclude = ppairs[(ppairs.event1 .== tpairs.event1[tn1maxi]) .&& (ppairs.event2 .== tpairs.event2[tn1maxi]), :]
  else
    ppairsexclude = ppairs[(ppairs.event1 .== tpairs.event1[tn2maxi]) .&& (ppairs.event2 .== tpairs.event2[tn2maxi]), :]
  end
else
  ppairsexclude = ppairs[(ppairs.event1 .== ppairsu.event1[pn2maxi]) .&& (ppairs.event2 .== ppairsu.event2[pn2maxi]), :]
end
ppairsexclude = ppairsexclude[:, [:station, :event1, :event2]]
rename!(ppairsexclude,:station => :pstation)
#excludepairs = vcat(excludepairs,ppairsexclude)
#global pn2max = ppairs.n2[pn2maxi]
#CSV.write("data/catalogs/japan_WAKE_exclude.csv", excludepairs)
tpairs.Δτ = vec(vcat([(tpairs.Δτ[i])' for i = 1:nt]...)) 
CSV.write(@sprintf("results/pairs/japan_%s_ppairsexclude.csv",tsname), ppairsexclude)
CSV.write(@sprintf("results/pairs/japan_%s_ppairs_2.5a3.5hz.csv",tsname), ppairsu)
CSV.write(@sprintf("results/pairs/japan_%s_tpairs_2.5a3.5hz.csv",tsname), tpairs)

hist = true
if hist
    b = -0.3:5e-3:0.3
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
#end

# reconstruct full travel time anomalies
τ = reshape(D*a, (m, 1)) 
τ2 = reshape(D*(E\y), (m, 1))
eτ = reshape(sqrt.(diag(D*P*D')), (m, 1))

filename = @sprintf("results/anomalies/japan_%s_x2y_40ex.h5",tsname)
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
CSV.write(@sprintf("results/anomalies/japan_%s.csv",tsname), df_T)

td = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24

# interpolate onto regular grid
θi = -9:0.25:9 
θd = θ

# interpolate onto regular grid
plot = false
for year = 1997:2021
  @printf("\ninterpolation year %d\n",year)
  ti = DateTime(year, 1, 1, 12, 0, 0) : Day(10) : DateTime(year, 12, 31, 12, 0, 0)
  τi,ei,ei0 = SOT.regulargrid(td, θd, ti, θi, a, R, invRyy, E, λt, λθ, στ)
  if year == 1997
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
#eiall = eiall./στ 
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
    
plot = false
if plot
    
    fig, ax = subplots(2,1,figsize=(190/25.4, 190/25.4),sharex=true)
    ax[1].pcolormesh(ti,θi,τi,vmin=-1,vmax=1,cmap="RdBu",shading="nearest",zorder=0)
    ax[1].pcolormesh(ti,θi,ei,vmin=0,vmax=1.5,cmap="Greys",alpha=0.2,shading="nearest")
    ax[2].pcolormesh(ti,θi,ri,vmin=0,vmax=1,shading="nearest")
    for i = 1:2
      ax[i].scatter(df_T.time,df_T.azimuth,s=5,color="black")
      ax[i].set_ylabel("azimuth (deg)")
      ax[i].set_xlim([ti[1],ti[end]])
    end
    fig.tight_layout()
    fig.savefig(@sprintf("results/plots/japan_%s_regularmap_2011.pdf",tsname))
end
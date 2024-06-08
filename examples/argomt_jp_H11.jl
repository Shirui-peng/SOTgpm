include("../src/SOT.jl")
using .SOT, Printf, LinearAlgebra, Dates, HDF5, DataFrames, CSV, Statistics
using Random, Distributions, SparseArrays, PyPlot, Interpolations, Cubature

# receiver data files
station = "H11N3"
nexclude = 16
ppfile = @sprintf("results/pairs/japan_%s_ppairs_2.5a3.5hz_%dex.csv",station,nexclude)
tpvfile = @sprintf("results/pairs/japan_%s_tpairs_2.5a3.5hz_%dex.csv",station,nexclude)

# receiver location
tstation = [166.90986,19.71786]

# collect and clean up pairs
tpairs = CSV.read(tpvfile, DataFrame)
ppairs = CSV.read(ppfile, DataFrame)
ppairs = innerjoin(ppairs, select(tpairs, [:event1, :event2]), on=[:event1, :event2])
ppairs.stn6 = [s[1:6] for s in ppairs.station]
unique!(ppairs, [:stn6,:event1,:event2])
select!(ppairs, Not([:stn6,:cc,:n1,:n2]))

# reference source location
evtpos = [142.85,38.10]

# p-wave station locations
pstations = ["IU.MAJO.00.BHZ","IU.MAJO.10.BHZ","PS.TSK..BHZ", "II.ERM.00.BHZ", "G.INU.00.BHZ"]
pstnlats,pstnlons = [36.55,36.55,36.21,42.02,35.35],[138.2,138.2,140.11,143.16,137.03]
pstations = DataFrame(station=pstations,slat=pstnlats,slon=pstnlons)

# get mean time
t1, = SOT.getE(tpairs,ppairs,tstation,pstations;hydro=true)
trd = Dates.value.(t1 - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
tm = (max(trd...)+min(trd...))/2

# path travel time anomaly trend prior (s/yr)
σtrend = 0.037

# path travel time anomaly annual cycle prior (s)
σannual = 0.1

# path travel time anomaly semi-annual cycle prior (s)
σsemiannual = 0.1

# time scale (d)
λt = [60,62,57]

# length scale (km)
λx = [100/sqrt(2),65,59]

# stochastic slowness anomaly scale (ms/km)
σs = [0.44,0.42,0.32]

# offset slowness anomaly scale (ms/km)
σm = 0.44

# noise (ms/km)
σn = [0,0.041,0.016]

# noise parameters for seismic data
σp,σx,σ,σnp,σh = 0.97,0.019,8.1e-3,2.2e-3,0.029

# stacked parameter vector
x = [λt[1],λx[1],σs[1],σm,σn[1],σp,σx,σ,σnp,σh]

#lone = h5read("results/anomalies/japan_H11N3_x2y_16ex.h5", "lon")
#late = h5read("results/anomalies/japan_H11N3_x2y_16ex.h5", "lat")
#azme = h5read("results/anomalies/japan_H11N3_x2y_16ex.h5", "θ")

# read conventional in situ data
ctdfile = "results/ship/wod/ctds_kuroshio.csv"
argofile = "results/argo/japan_1997-01_2022-12.csv"
floats = vcat(CSV.read(ctdfile, DataFrame)[:,1:7],CSV.read(argofile, DataFrame)[:,1:7])

# constrian in situ data manually
t0,t1,y0,z0 = DateTime(2007, 12, 16),DateTime(2021, 8, 1),30,-1.9e3
floats = floats[(floats.y .<=36) .| (floats.x .>= 141),:]
floats = floats[(floats.z0.<z0),:]
floats = floats[(t0 .≤ floats.t .< t1),:]
floats = floats[.!(isnan.(floats.Δsp1)),:]
floats = floats[abs.(floats.Δsp1) .< 10,:]
floats = floats[floats.zmax.>-1.9e3,:]
unique!(floats, [:x,:y,:t])
mf = size(floats,1)
println(mf)

# regular-grid path-path covariance calculation
@printf("covariance integration...time: %s\n",Dates.format(now(), "HH:MM:SS"))
θgrid = -0.1:0.05:20
cgrid = SOT.cint.(θgrid)
itp_c = linear_interpolation(θgrid, cgrid, extrapolation_bc=0)
θgridt = -0.2:0.05:20
cgridt = SOT.tcint.(θgridt)
itp_tc = linear_interpolation(θgridt, cgridt, extrapolation_bc=0)
@printf("path-path covariance doen. time: %s\n",Dates.format(now(), "HH:MM:SS"))

# regular-grid point-path covariance calculation
θgrid2 = -10:0.1:10
xgrid,ygrid = 135:0.1:170,15:0.1:45
ctgrid = SOT.cinta.(xgrid,ygrid',reshape(θgrid2,1,1,:))
nodes = (xgrid, ygrid, θgrid2)
itp_ca = interpolate(nodes, ctgrid, Gridded(Linear()))
θgrid2t = -10:0.1:10
xgridt,ygridt = 135:0.1:170,15:0.1:45
nodest = (xgridt, ygridt, θgrid2t)
tctgrid = SOT.tcinta.(xgridt,ygridt',reshape(θgrid2t,1,1,:))
itp_tca = interpolate(nodest, tctgrid, Gridded(Linear()))
@printf("point-path covariance doen. time: %s\n",Dates.format(now(), "HH:MM:SS"))

# read map grid from rg product
h5fl = "results/argo/kuroshio_rg.h5"
xrg = h5read(h5fl, "x")
yrg = h5read(h5fl, "y")
Δsrg = h5read(h5fl, "Δsp1")[:,:,1]
xidx = 135 .< xrg .< 170
nxrg,nyrg=length(xrg),length(yrg)

# initialze offset and trend map
gmn,gtd = fill(NaN, nxrg, nyrg, 2),fill(NaN, nxrg,nyrg, 6)

# iterate over latitude 
for (j,y) in enumerate(yrg)
  if 15<y<45
  
    # select in situ data within 3-deg latitude 
    floatsj = floats[abs.(floats.y.-y).<3,:]
    @printf("lat %.2f, npoint %d\n",y,size(floatsj,1))
    
    # select grid points to estimate offset and trend
    idx = .!iszero.(Δsrg[:,j]) .& xidx
    if y>36
      idx = idx .& (xrg.>141)
    end
    nxj = sum(idx)
    gxy = [xrg[idx] y*ones(nxj)]
    
    # estimate offset and trend with seismic, conventional, and combined data
    Pa,xa,Pt,xt,Pat,xat = SOT.invertargotmt(gxy,x,σtrend,σannual,σsemiannual,tm,floatsj,tpairs,ppairs,tstation,pstations,evtpos,itp_c,itp_tc,itp_ca,itp_tca)
    gmn[idx,j,1],gtd[idx,j,1]=xa[1:nxj],xa[nxj+1:end]
    gmn[idx,j,2],gtd[idx,j,2]=sqrt.(diag(Pa)[1:nxj]),sqrt.(diag(Pa)[nxj+1:end])
    gtd[idx,j,3],gtd[idx,j,4]=xt[nxj+1:end],sqrt.(diag(Pt)[nxj+1:end])
    gtd[idx,j,5],gtd[idx,j,6]=xat[nxj+1:end],sqrt.(diag(Pat)[nxj+1:end])
  end
end

# save to file
h5open("results/argo/kuroshio_mt_H11.h5", "w") do file
  write(file, "x", xrg)
  write(file, "y", yrg)
  write(file, "gmn", gmn)
  write(file, "gtd", gtd)
end
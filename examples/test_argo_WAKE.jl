include("../src/SOT.jl")
using .SOT, Printf, LinearAlgebra, Dates, HDF5, DataFrames, CSV, Statistics
using Random, Distributions, SparseArrays, PyPlot, Interpolations, Cubature

# receiver location
tstation = [166.652,19.283]

# reference source location
evtpos = [142.85,38.10]

# reference azimuth
θ0 = SOT.azimuth(tstation[2],tstation[1],evtpos[2],evtpos[1])

# azimuth range
θ1,θ2=-6.5,3.5

# source location on path 
h11a0 = [142.06,38.87]

# events for inversion
te = DateTime(1997, 7, 1):Day(10):DateTime(2021, 7, 31)
m = length(te)
lone = ones(m)*h11a0[1]
late = ones(m)*h11a0[2]
azme = 0.6*ones(m)
nms = ["longitude","latitude","time","azimuth"]
events = DataFrame([lone late te azme], nms)

# receiver data files
station = "WAKE"
nexclude = 40
ppfile = @sprintf("results/pairs/japan_%s_ppairs_2.5a3.5hz_%dex.csv",station,nexclude)
tpvfile = @sprintf("results/pairs/japan_%s_tpairs_2.5a3.5hz_%dex.csv",station,nexclude)

# collect and clean up pairs
tpairs = CSV.read(tpvfile, DataFrame)
ppairs = CSV.read(ppfile, DataFrame)
tt = DateTime(2000, 1, 1) .+ Millisecond.(h5read("results/anomalies/japan_WAKE_x2y_40ex.h5", "t"))
θt = h5read("results/anomalies/japan_WAKE_x2y_40ex.h5", "θ")
ia = θ1.<=θt.<=θ2
tt = tt[ia]
#tpairs = filter(x -> (x[:event1] in tt) & (x[:event2] in tt), tpairs)
unique!(tpairs, [:event1,:event2])
ppairs = innerjoin(ppairs, select(tpairs, [:event1, :event2]), on=[:event1, :event2])
ppairs.stn6 = [s[1:6] for s in ppairs.station]
unique!(ppairs, [:stn6,:event1,:event2])
select!(ppairs, Not([:stn6,:cc,:n1,:n2]))

# p-wave station locations
pstations = ["IU.MAJO.00.BHZ","IU.MAJO.10.BHZ","PS.TSK..BHZ", "II.ERM.00.BHZ", "G.INU.00.BHZ"]
pstnlats,pstnlons = [36.55,36.55,36.21,42.02,35.35],[138.2,138.2,140.11,143.16,137.03]
pstations = DataFrame(station=pstations,slat=pstnlats,slon=pstnlons)

# pair sizes
nt1 = size(tpairs, 1)
np1 = size(ppairs, 1)

# path travel time anomaly trend prior (s/yr)
σtrend = 0.037

# path travel time anomaly annual cycle prior (s)
σannual = 0.1

# path travel time anomaly semi-annual cycle prior (s)
σsemiannual = 0.1

# time scale (d)
λt = [60,64,57]

# length scale (km)
λx = [100/sqrt(2),65,59]

# stochastic slowness anomaly scale (ms/km)
σs = [0.44,0.50,0.32]

# offset slowness anomaly scale (ms/km)
σm = 0.44

# noise (ms/km)
σn = [0.0,0.049,0.016]

# noise parameters for seismic data
σp,σx,σ,σnp,σh = 0.86,0.028,10e-3,2.9e-3,0

# stacked parameter vector
i = 1
x = [λt[i],λx[i],σs[i],σm,σn[i],σp,σx,σ,σnp,σh]

# regular-grid path-path covariance calculation
@printf("covariance integration...time: %s\n",Dates.format(now(), "HH:MM:SS"))
θgrid = -0.1:0.05:20
cgrid = SOT.cint.(θgrid)
itp_c = linear_interpolation(θgrid, cgrid, extrapolation_bc=0)
θgridt = -0.2:0.05:20
cgridt = SOT.tcint.(θgridt)
itp_tc = linear_interpolation(θgridt, cgridt, extrapolation_bc=0)
@printf("path covariance done. time: %s\n",Dates.format(now(), "HH:MM:SS"))

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
@printf("point-path covariance done. time: %s\n",Dates.format(now(), "HH:MM:SS"))

# read conventional in situ data
ctdfile = "results/ship/wod/ctds_kuroshio.csv"
argofile = "results/argo/japan_1997-01_2022-12.csv"
floats = vcat(CSV.read(ctdfile, DataFrame)[:,1:7],CSV.read(argofile, DataFrame)[:,1:7])

# constrian in situ data manually
t0,t1,x0,y0,z0 = DateTime(1997,6,5),DateTime(2021,8,1),141,36,-1.9e3
floats = floats[(floats.z0.<z0),:]
floats = floats[(floats.y .<=y0) .| (floats.x .>= x0),:]
floats = floats[(t0 .≤ floats.t .< t1),:]
floats = floats[.!(isnan.(floats.Δsp1)),:]
floats = floats[abs.(floats.Δsp1) .< 10,:]
floats = floats[floats.zmax.>z0,:]
unique!(floats, [:x,:y,:t])
d0 = 300
#Δt = 183
nxk = 321
#nearfloats = SOT.nearargo(floats,events,nxk,d0,Δt,tstation)
xk, yk = SOT.findpath(tstation, h11a0, nxk)
idxk = minimum(SOT.dist.(xk', yk', floats.x, floats.y),dims=2)./1e3 .< d0
nearfloats = floats[vec(idxk),:]
#idxa = θ1.<= SOT.azimuth.(tstation[2],tstation[1],floats.y, floats.x) .- θ0 .<= θ2
#nearfloats = floats[vec(idxa),:]
println(size(nearfloats,1))

# invert for path-mean temperatures
Pa,xa,Pt,xt,Pat,xat,R0,dR = SOT.invertargot0(te,azme,x,σtrend, σannual, σsemiannual,tpairs,ppairs,tstation,evtpos,pstations,nearfloats,itp_c,itp_ca,itp_tc,itp_tca)
#SOT.invertargottc(te,azme,x,3.7σtrend, σannual, σsemiannual,tpairs,ppairs,tstation,evtpos,pstations,nearfloats,λx[1],σs[1],σn[1],itp_c,itp_ca,itp_tc,itp_tca)

# save to file
m = length(te)
ω = 2π/SOT.meanyear
ter = Dates.value.(te - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
tm = (max(ter...)+min(ter...))/2
D = [I(m) ter.-tm cos.(ω*ter) sin.(ω*ter) cos.(2ω*ter) sin.(2ω*ter)]
h5open("results/argo/point_WAKE_.6d.h5", "w") do file
  write(file, "xa", xa)
  write(file, "xt", xt)
  write(file, "xat", xat)
  write(file, "Pa", Pa)
  write(file, "Pt", Pt)
  write(file, "Pat", Pat)
  write(file, "e0", sqrt.(diag(Matrix(D*R0*D'))))
  write(file, "t", Dates.value.(te - DateTime(2000, 1, 1)))
end

include("../src/SOT.jl")
using .SOT, Printf, LinearAlgebra, Dates, HDF5, DataFrames, CSV, Statistics
using Random, Distributions, SparseArrays, PyPlot, Interpolations, HCubature

# event times for inversion
te = DateTime(1997, 7, 1):Day(10):DateTime(2021, 7, 31)
m = length(te)

# receiver location
tstation = [166.652,19.283]

# reference source location
evtpos = [142.85,38.10]

# reference azimuth
θ0 = SOT.azimuth(tstation[2],tstation[1],evtpos[2],evtpos[1])

# area shape
shape = "tri"

# azimuth and length limits
θ1,θ2,L=-6.5,3.5,3.2

# read seismic event data
h5fl = "results/anomalies/japan_WAKE_x2y_40ex.h5"
lone = h5read(h5fl, "lon")
late = h5read(h5fl, "lat")
azme = h5read(h5fl, "θ")
tevt = DateTime(2000, 1, 1) .+ Millisecond.(h5read(h5fl, "t"))
tevt = tevt[θ1.<=azme.<=θ2]

# receiver data files
station = "WAKE"
nexclude = 40
ppfile = @sprintf("results/pairs/japan_%s_ppairs_2.5a3.5hz_%dex.csv",station,nexclude)
tpvfile = @sprintf("results/pairs/japan_%s_tpairs_2.5a3.5hz_%dex.csv",station,nexclude)

# collect and clean up pairs
tpairs = CSV.read(tpvfile, DataFrame)
ppairs = CSV.read(ppfile, DataFrame)
tpairs = filter(x -> (x[:event1] in tevt) & (x[:event2] in tevt), tpairs)
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
nxk = 321
#is5,in5 = argmin(azme.+1.5),argmin(azme.-5)
#xks, yks = SOT.findpath(tstation, (lone[is5],late[is5]), nxk)
#xkn, ykn = SOT.findpath(tstation, (lone[in5],late[in5]), nxk)
#d0 = 300
#idxs = minimum(SOT.dist.(xks', yks', floats.x, floats.y),dims=2)./1e3 .< d0
#idxn = minimum(SOT.dist.(xkn', ykn', floats.x, floats.y),dims=2)./1e3 .< d0
idxa = θ1.<= SOT.azimuth.(tstation[2],tstation[1],floats.y, floats.x) .- θ0 .<= θ2
floats = floats[vec(idxa),:]# .| vec(idxs),:]
println(size(floats,1))
# real time (days)
mf = size(floats,1)
#floats.Δsp1 .-= mean(floats.Δsp1)

println("covariance integration...")
itp = true
if itp

    # regular-grid path-path covariance calculation
    θgrid = -1:0.05:20
    cgrid = SOT.cint.(θgrid)
    itp_c = linear_interpolation(θgrid, cgrid, extrapolation_bc=0)
    cgrid = SOT.tcint.(θgrid)
    itp_tc = linear_interpolation(θgrid, cgrid, extrapolation_bc=0)
    
    # regular-grid point-path covariance calculation
    θgrid = -10:0.1:10
    xgrid,ygrid = 135:0.1:170,15:0.1:45
    ctgrid = SOT.cinta.(xgrid,ygrid',reshape(θgrid,1,1,:))
    nodes = (xgrid, ygrid, θgrid)
    itp_ca = interpolate(nodes, ctgrid, Gridded(Linear()))
    tctgrid = SOT.tcinta.(xgrid,ygrid',reshape(θgrid,1,1,:))
    nodes = (xgrid, ygrid, θgrid)
    itp_tca = interpolate(nodes, tctgrid, Gridded(Linear()))
    
    # regular-grid path-area covariance calculation
    θgrid = -10:0.05:10
    cgrid = SOT.cint2d.(θgrid;shape)
    itp_c2d = linear_interpolation(θgrid, cgrid, extrapolation_bc=0)
    cgrid = SOT.tcint2d.(θgrid;shape)
    itp_tc2d = linear_interpolation(θgrid, cgrid, extrapolation_bc=0)
    
     # regular-grid point-area covariance calculation
    nodes = (xgrid, ygrid)
    c2grid = SOT.cint2da.(xgrid,ygrid';shape)
    itp_ca2d = interpolate(nodes, c2grid, Gridded(Linear()))
    c2grid = SOT.tcint2da.(xgrid,ygrid';shape)
    itp_tca2d = interpolate(nodes, c2grid, Gridded(Linear()))
end

# area-mean auto covariance calculation
if shape=="tri"
    c0(x) = sqrt(2)*x[1]*x[2]*exp(-sqrt(x[1]^2+x[2]^2-2x[1]*x[2]*cosd(x[3]-x[4]))/(1e-3λx[1])/sqrt(2))*cos(sqrt(x[1]^2+x[2]^2-2x[1]*x[2]*cosd(x[3]-x[4]))/(1e-3λx[1])/sqrt(2)-π/4)
    c2d = hcubature(c0, [0,0,θ1,θ1], [L,L,θ2,θ2]; rtol=1e-6, atol=1e-6)[1]/(((θ2-θ1)*L/2)^2)
    c0(x) = x[1]*x[2]*exp.(-sqrt(x[1]^2+x[2]^2-2x[1]*x[2]*cosd(x[3]-x[4]))/0.2)
    tc2d = hcubature(c0, [0,0,θ1,θ1], [L,L,θ2,θ2]; rtol=1e-6, atol=1e-6)[1]/(((θ2-θ1)*L/2)^2)
else
    c0(x) = sqrt(2)*exp(-sqrt((x[1]-x[2])^2+(x[3]-x[4])^2)/(1e-3λx[1])/sqrt(2))*cos(sqrt((x[1]-x[2])^2+(x[3]-x[4])^2)/(1e-3λx[1])/sqrt(2)-π/4)
    c2d = hcubature(c0, [0,0,π*L*θ1/180,π*L*θ1/180], [L,L,π*L*θ2/180,π*L*θ2/180]; rtol=1e-6, atol=1e-6)[1]/((π*L*(θ2-θ1)/180)^2)
    c0(x) = exp.(-sqrt((x[1]-x[2])^2+(x[3]-x[4])^2)/0.2)
    tc2d = hcubature(c0, [0,0,π*L*θ1/180,π*L*θ1/180], [L,L,π*L*θ2/180,π*L*θ2/180]; rtol=1e-6, atol=1e-6)[1]/((π*L*(θ2-θ1)/180)^2)
end
println("covariance integration finished.")

# invert for area-mean temperature anomalies
Pa,xa,Pt,xt,Pat,xat,C0 = SOT.invertargot2dtc(te,x,3.7σtrend, σannual, σsemiannual,tpairs,ppairs,tstation,evtpos,pstations,floats,c2d,tc2d;
                                             itp_c,itp_ca,itp_tc,itp_tca,itp_c2d,itp_ca2d,itp_tc2d,itp_tca2d,full=true)

# save to file
ter = Dates.value.(te - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
tm = (max(ter...)+min(ter...))/2
ω = 2π/SOT.meanyear
m = length(te)
D = [I(m) ter.-tm cos.(ω*ter) sin.(ω*ter) cos.(2ω*ter) sin.(2ω*ter)]
h5open(@sprintf("results/argo/point_WAKE_2d%s.h5",shape), "w") do file
  write(file, "xa", xa)
  write(file, "xt", xt)
  write(file, "xat", xat)
  write(file, "Pa", Pa)
  write(file, "Pt", Pt)
  write(file, "Pat", Pat)
  write(file, "e0", sqrt.(diag(Matrix(D*C0*D'))))
  write(file, "c0", sqrt.(diag(Matrix(C0))))
  write(file, "t", Dates.value.(te - DateTime(2000, 1, 1)))
end

include("../src/SOT.jl")
using PyCall, PyPlot, Printf, Dates, LinearAlgebra, Statistics
using SparseArrays, DataFrames, CSV, Trapz
using HDF5, NCDatasets, Interpolations, Random, Distributions

# Earth's radius
const earthradius = 6371e3

gsw = pyimport("gsw")
np = pyimport("numpy")
argo = pyimport("argopy")
loader = argo.DataFetcher()

"Calculate great-circle distance"
function dist(lon1, lat1, lon2, lat2)
  cosΔσ = sind(lat1)*sind(lat2) + cosd(lat1)*cosd(lat2)*cosd(lon1-lon2)
  return earthradius*acos(sign(cosΔσ)*min(1,abs(cosΔσ)))
end

"in situ temperature from potential temperature" 
function T(Θ, S, λ, θ, z)
  if ismissing(Θ)
    return missing
  else
    p = gsw.p_from_z(z, θ)
    SA = gsw.SA_from_SP(S, p, λ, θ)
    CT = gsw.CT_from_pt(SA, Θ)
    T = gsw.t_from_CT(SA, CT, p)
    return T
  end
end

"""
lon. and lat. of n equally spaced points along the great circle path that 
originates at p1d = (λ1, θ1) and passes through p2d = (λ2, θ2)
"""
function findpath(p1d, p2d, n)
  λ1 = deg2rad(p1d[1])
  θ1 = deg2rad(p1d[2])
  λ2 = deg2rad(p2d[1])
  θ2 = deg2rad(p2d[2])
  λ2p = atan(cos(θ2)*sin(λ2-λ1), sin(θ1)*cos(θ2)*cos(λ2-λ1) - cos(θ1)*sin(θ2))
  θ2p = asin(cos(θ1)*cos(θ2)*cos(λ2-λ1) + sin(θ1)*sin(θ2))
  R2 = [cos(λ2p) -sin(λ2p) 0; sin(λ2p) cos(λ2p) 0; 0 0 1]
  Ry = [sin(θ1) 0 cos(θ1); 0 1 0; -cos(θ1) 0 sin(θ1)]
  Rz = [cos(λ1) -sin(λ1) 0; sin(λ1) cos(λ1) 0; 0 0 1]
  θp = π/2 .- (0:n-1)/(n-1)*(π/2 - θ2p)
  p = [cos.(θp)'; zeros(n)'; sin.(θp)']
  q = Rz*Ry*R2*p
  λq = atan.(q[2,:], q[1,:])
  θq = asin.(q[3,:])
  return rad2deg.(λq), rad2deg.(θq)
end

# load coordinates and variables
xe = h5read("results/ecco/wp_kuroshio_0817.h5", "x")
ye = h5read("results/ecco/wp_kuroshio_0817.h5", "y")
ze = h5read("results/ecco/wp_kuroshio_0817.h5", "z")
wp = h5read("results/ecco/wp_kuroshio_0817.h5", "wp")
H = h5read("results/ecco/wp_kuroshio_0817.h5", "H")
#G = h5read("results/ecco/G_kuroshio_0817.h5", "G")
#n2 = h5read("results/ecco/G_kuroshio_0817.h5", "n2")

fl = Dataset("../argo/temp/RG_ArgoClim_Temperature_2019.nc", "r")
xa = Array{Float64}(fl["LONGITUDE"][:])
ya = Array{Float64}(fl["LATITUDE"][:])
pm = Array{Float64}(fl["PRESSURE"][:])

xmin,xmax,ymin,ymax = 135,170,15.,45.
xidx = findall(xmin .<= xa .<= xmax)
yidx = findall(ymin .<= ya .<= ymax)
#Ta = Array{Float64}(replace(fl["ARGO_TEMPERATURE_MEAN"][:,:,:], missing=>NaN))
#ΔTa = Array{Float64}(replace(fl["ARGO_TEMPERATURE_ANOMALY"][xidx,yidx,:,:], missing=>0))
xa,ya = xa[xidx],ya[yidx]
nxa,nya = length(xa),length(ya)
#ta = DateTime(2004, 1, 1):Month(1):DateTime(2019, 1, 1)
ta = DateTime(2019, 1, 1):Month(1):DateTime(2022, 1, 1)
ta = ta[1:end-1] .+ 0.5*(ta[2:end] .- ta[1:end-1])
nta = length(ta)

na = length(pm)
# node locations (made-up non-uniform grid)
pa = zeros(na+1)

for k = 1:na
  pa[1+k] = 2*pm[k]-pa[k]
end

iR = inv(Diagonal([1, 0.5, 0.25, 0.1].^2))
iN = 0.1^-2*I

# mass matrix
Im = Int[]; J = Int[]; V = Float64[]
for i = 1:na
  push!(Im, i); push!(J, i); push!(V, (pa[i+1]-pa[i])/3)
  push!(Im, i); push!(J, i+1); push!(V, (pa[i+1]-pa[i])/6)
  push!(Im, i+1); push!(J, i); push!(V, (pa[i+1]-pa[i])/6)
  push!(Im, i+1); push!(J, i+1); push!(V, (pa[i+1]-pa[i])/3)
end
M = sparse(Im, J, V)

# conversion matrix
Ic = Int[]; J = Int[]; V = Float64[]
for i = 1:na
push!(Ic, i); push!(J, i); push!(V, (pa[i+1]-pa[i])/2)
push!(Ic, i+1); push!(J, i); push!(V, (pa[i+1]-pa[i])/2)
end
C = sparse(Ic, J, V)

nmode = 4
Δsp1 = zeros(nxa,nya,nta)
Δsp2 = zeros(nxa,nya,nta)
for n = 1:nta
    @printf("projection for %s\n",ta[n])
    flnm = @sprintf("../argo/sio-argo.ucsd.edu/pub/www-argo/RG/RG_ArgoClim_%4d%02d_2019.nc",Dates.year(ta[n]),Dates.month(ta[n]))
    nfl = Dataset(flnm, "r")
    ΔTa = Array{Float64}(replace(Dataset(flnm, "r")["ARGO_TEMPERATURE_ANOMALY"][xidx,yidx,:,:], missing=>0))
    for j = 1:nya
      je = argmin(abs.(ye.-ya[j]))
      z = gsw.z_from_p.(pa, ya[j])
      #@printf("projection for latitude %.1f\n",ya[j])
      for i = 1:nxa
        ie = argmin(abs.(xe.-xa[i]))
        if H[ie,je]>2e3
          ΔT = M\(C*ΔTa[i,j,:,1])
          etpw1 = linear_interpolation(ze, wp[ie,je,:,1],extrapolation_bc=Line()) 
          #Δsp1[i,j,:] = -1e6*[trapz(z,ΔT[:,n].*etpw1.(z)) for n = 1:nta]
          Δsp1[i,j,n] = -1e6*trapz(z,ΔT.*etpw1.(z))
          #etpw2 = linear_interpolation(ze, wp[ie,je,:,2],extrapolation_bc=Line()) 
          #Δsp2[i,j,:] = -1e6*[trapz(z,ΔT[:,n].*etpw2.(z)) for n = 1:nta]
        end
      end
    end
end

jy,ix = argmin(abs.(ya.-35)),argmin(abs.(xa.-145))
fig,ax=subplots(figsize=(7.2,3.2))
ax.plot(ta,Δsp1[ix,jy,:])
fig.savefig("results/argo/kuroshio_rgtest.pdf")

h5open("results/argo/kuroshio_rg1922.h5", "w") do file
  write(file, "x", xa)
  write(file, "y", ya)
  write(file, "t", Dates.value.(ta - DateTime(2000, 1, 1)))
  write(file, "Δsp1", Δsp1)
  write(file, "Δsp2", Δsp2)
end
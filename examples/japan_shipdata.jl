include("../src/SOT.jl")
using PyCall, PyPlot, Printf, Dates, LinearAlgebra, Statistics
using SparseArrays, DataFrames, CSV, Trapz, Glob
using HDF5, NCDatasets, Interpolations, Random, Distributions

# Earth's radius
const earthradius = 6371e3

gsw = pyimport("gsw")
np = pyimport("numpy")

"Calculate great-circle distance"
function dist(lon1, lat1, lon2, lat2)
  cosΔσ = sind(lat1)*sind(lat2) + cosd(lat1)*cosd(lat2)*cosd(lon1-lon2)
  return earthradius*acos(sign(cosΔσ)*min(1,abs(cosΔσ)))
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

function getMC(da)
  na = length(da)-1
  # mass matrix
  Im = Int[]; J = Int[]; V = Float64[]
  for i = 1:na
    push!(Im, i); push!(J, i); push!(V, (da[i+1]-da[i])/3)
    push!(Im, i); push!(J, i+1); push!(V, (da[i+1]-da[i])/6)
    push!(Im, i+1); push!(J, i); push!(V, (da[i+1]-da[i])/6)
    push!(Im, i+1); push!(J, i+1); push!(V, (da[i+1]-da[i])/3)
  end
  M = sparse(Im, J, V)

  # conversion matrix
  Ic = Int[]; J = Int[]; V = Float64[]
  for i = 1:na
    push!(Ic, i); push!(J, i); push!(V, (da[i+1]-da[i])/2)
    push!(Ic, i+1); push!(J, i); push!(V, (da[i+1]-da[i])/2)
  end
  C = sparse(Ic, J, V)
  return M,C
end

nxk,d0 = 321,500
tstation = (166.90986,19.71786)
evtpos = (142.85,38.10)
θ0 = SOT.azimuth(tstation[2],tstation[1],evtpos[2],evtpos[1])
xt = h5read("results/anomalies/japan_H11N3_x2y_16ex.h5", "lon")
yt = h5read("results/anomalies/japan_H11N3_x2y_16ex.h5", "lat")
at = h5read("results/anomalies/japan_H11N3_x2y_16ex.h5", "θ")
iamin,iamax = argmin(at),argmax(at)
xks, yks = findpath(tstation, (xt[iamin],yt[iamin]), nxk)
xkn, ykn = findpath(tstation, (xt[iamax],yt[iamax]), nxk)

# load coordinates and variables
xe = h5read("results/ecco/wp_kuroshio_0817.h5", "x")
ye = h5read("results/ecco/wp_kuroshio_0817.h5", "y")
ze = h5read("results/ecco/wp_kuroshio_0817.h5", "z")
wp = h5read("results/ecco/wp_kuroshio_0817.h5", "wp")
H = h5read("results/ecco/wp_kuroshio_0817.h5", "H")
Te = h5read("results/ecco/wp_kuroshio_0817.h5", "T")
G = h5read("results/ecco/G_kuroshio_0817.h5", "G")
n2 = h5read("results/ecco/G_kuroshio_0817.h5", "n2")

woaclim = true

if woaclim
  xmin,xmax,ymin,ymax=130,170,10,50
  fl = Dataset("data/ship/woa18_A5B7_t00_01.nc", "r")
  xa = Array{Float64}(fl["lon"][:])
  ya = Array{Float64}(fl["lat"][:])
  da = Array{Float64}(fl["depth_bnds"][:,:])
  da = [da[1,:]; da[2,end]]
  xidx = xmin .< xa .< xmax
  yidx = ymin .< ya .< ymax
  Ta = Array{Float64}(replace(fl["t_mn"][xidx,yidx,:,1], missing=>NaN))
  xa,ya = xa[xidx],ya[yidx]
end

# list of all selected T-wave pairs
ctds = DataFrame(x=Float64[], y=Float64[], z0=Float64[], t=DateTime[],
                 Δsp1=Float64[], Δsp2=Float64[], zmax=Float64[], smax=Float64[],
                 Δsp1r=Float64[], Δsp2r=Float64[])

filenames = glob("wod_*.nc","data/ship/wod")
plot = false
#nf = length(filenames)
rc("font", size=8)
rc("axes", titlesize="medium")
for fn in filenames
  fl = Dataset(fn, "r")
  x = [Float64(fl["lon"][:])]#Array{Float64}(fl["longitude"][:])
  y = [Float64(fl["lat"][:])]#Array{Float64}(fl["latitude"][:])
  t = [DateTime(fl["time"][:])]#Array{DateTime}(fl["time"][:])
  z = -Array{Float64}(replace(fl["z"][:], missing=>NaN))
  z_qc = Array{Int16}(replace(fl["z_WODflag"][:], missing=>9))
  #p = Array{Float64}(replace(fl["pressure"][:,:], missing=>NaN))
  #p_qc = Array{Int8}(replace(fl["pressure_qc"][:,:], missing=>9))
  T = Array{Float64}(fl["Temperature"][:])#Array{Float64}(replace(fl["ctd_temperature"][:,:], missing=>NaN))
  T_qc = Array{Int16}(fl["Temperature_WODflag"][:])#Array{Int8}(replace(fl["ctd_temperature_qc"][:,:], missing=>9))
  #depth = Array{Float64}(replace(fl["btm_depth"][:], missing=>NaN))
  J = 1#length(x)
  idxk = (z_qc.==0) .& (T_qc.==0)
  for j = 1:J
    if sum(idxk)>0#!isnan(depth[j])
      ie,je = argmin(abs.(xe.-x[j])),argmin(abs.(ye.-y[j]))
      #idxk = (p_qc[:,j].==2) .& (T_qc[:,j].==2)
      #z = gsw.z_from_p.(p[idxk,j], y[j])
      z = z[idxk]
      try
        ik0 = argmin(z)
        if z[ik0]<-2e3 && H[ie,je]>2e3
          #dists = minimum(SOT.dist.(xks, yks, x[j], y[j]))/1e3
          #distn = minimum(SOT.dist.(xkn, ykn, x[j], y[j]))/1e3
          #azm = SOT.azimuth(tstation[2],tstation[1],y[j],x[j])-θ0
          if true#at[iamin] < azm < at[iamax] || dists<d0 || distn<d0
            ia,ja = argmin(abs.(xa.-x[j])),argmin(abs.(ya.-y[j]))
            Taij = Ta[ia,ja,:]
            ka = findlast(.!isnan.(Taij))
            M,C = getMC(da[1:ka+1])
            etpT = linear_interpolation(da[1:ka+1], M\(C*Taij[1:ka]),extrapolation_bc=Line()) 
            ΔT = T[idxk] .- etpT.(-z)
            ΔT[z.<-da[ka+1]] .= 0
            etpw1 = linear_interpolation(ze, wp[ie,je,:,1],extrapolation_bc=Line()) 
            etpw2 = linear_interpolation(ze, wp[ie,je,:,2],extrapolation_bc=Line()) 
            Δsp1 = ΔT.*etpw1.(z)
            kmax = argmax(abs.(Δsp1))
            smax = 1e9*Δsp1[kmax]
            Δsp1 = -1e6*trapz(z,ΔT.*etpw1.(z))
            Δsp2 = -1e6*trapz(z,ΔT.*etpw2.(z))
            klow = z .< -2e3
            Δsp1r = -1e6*trapz(z[klow],ΔT[klow].*etpw1.(z[klow]))
            Δsp2r = -1e6*trapz(z[klow],ΔT[klow].*etpw2.(z[klow]))

            if !isnan(Δsp1)
              push!(ctds,[x[j],y[j],z[ik0],t[j],Δsp1,Δsp2,z[kmax],smax,Δsp1r,Δsp2r])
              #@printf("%d near ctds\n",size(ctds,1))
              CSV.write(@sprintf("results/ship/wod/ctds_%s.csv","kuroshio"), ctds)
              
              #if mod1(size(ctds,1),300)==3 && plot
              if Δsp1<-2.2 && plot
                @printf("%s has outlying slowness anomaly!\n",fn)
                fig,ax = subplots(1,3,sharey=true)
                ax[1].plot(ΔT,z/1e3,label="data")
                ax[1].set_xlabel("temperature anomaly (K)")
                ax[1].set_ylabel("z (km)")
                ax[2].plot(1e9*wp[ie,je,:,1],ze/1e3,label = "mode 1")
                ax[2].plot(1e9*wp[ie,je,:,2],ze/1e3,label = "mode 2")
                ax[2].legend(frameon=false)
                ax[2].set_xlabel("mode (s\$~\\mathrm{km}^{-2}~\\mathrm{K}^{-1}\$)")
                ax[3].plot(1e9*ΔT.*etpw1.(z),z/1e3)
                ax[3].plot(1e9*ΔT.*etpw2.(z),z/1e3)
                ax[3].set_xlabel("projection (s\$~\\mathrm{km}^{-2}\$)")
                fig.suptitle(@sprintf("x=%.1f, y=%.1f, t=%s, Δsp1=%.2f",x[j],y[j],t[j],Δsp1))
                fig.tight_layout()
                fig.savefig(@sprintf("results/ship/wod/ctd%.0fx%.0fy%s_zT%s.pdf",x[j],y[j],t[j],fn[end-13:end-3]))
                close(fig)
              end
            end
          end
        end
      catch e
        @printf("Error at x=%.1f, y=%.1f, t=%s\n",x[j],y[j],t[j])
        println(e)
      end
    end
  end
end

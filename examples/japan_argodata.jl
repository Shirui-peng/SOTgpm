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

nxk = 321
tstation = (166.90986,19.71786)
evtpos = (142.85,38.10)
θ0 = SOT.azimuth(tstation[2],tstation[1],evtpos[2],evtpos[1])
xt = h5read("results/pairs/japan_H11N3_x2y.h5", "lon")
yt = h5read("results/pairs/japan_H11N3_x2y.h5", "lat")
at = h5read("results/pairs/japan_H11N3_x2y.h5", "θ")
iamin,iamax = argmin(at),argmax(at)
xks, yks = findpath(tstation, (xt[iamin],yt[iamin]), nxk)
xkn, ykn = findpath(tstation, (xt[iamax],yt[iamax]), nxk)

# load coordinates and variables
xe = h5read("results/ecco/wp_kuroshio_0817.h5", "x")
ye = h5read("results/ecco/wp_kuroshio_0817.h5", "y")
ze = h5read("results/ecco/wp_kuroshio_0817.h5", "z")
wp = h5read("results/ecco/wp_kuroshio_0817.h5", "wp")
H = h5read("results/ecco/wp_kuroshio_0817.h5", "H")
G = h5read("results/ecco/G_kuroshio_0817.h5", "G")
n2 = h5read("results/ecco/G_kuroshio_0817.h5", "n2")

rgclim = false

if rgclim
  fl = Dataset("../argo/temp/RG_ArgoClim_Temperature_2019.nc", "r")
  xa = Array{Float64}(fl["LONGITUDE"][:])
  ya = Array{Float64}(fl["LATITUDE"][:])
  pm = Array{Float64}(fl["PRESSURE"][:])
  Ta = Array{Float64}(replace(fl["ARGO_TEMPERATURE_MEAN"][:,:,:], missing=>NaN))
end

woaclim = true

if woaclim
  xmin,xmax,ymin,ymax=130,170,10,50
  fl = Dataset("data/ship/woa18_A5B7_t00_01.nc", "r")
  xa = Array{Float64}(fl["lon"][:])
  ya = Array{Float64}(fl["lat"][:])
  dpa = Array{Float64}(fl["depth_bnds"][:,:])
  dpa = [dpa[1,:]; dpa[2,end]]
  xidx = xmin .< xa .< xmax
  yidx = ymin .< ya .< ymax
  Ta = Array{Float64}(replace(fl["t_mn"][xidx,yidx,:,1], missing=>NaN))
  xa,ya = xa[xidx],ya[yidx]
end

time_windows = [["1997-01-01","1999-12-31"]]
quarters = [["01-01", "03-31"], ["04-01", "06-30"], ["07-01", "09-30"], ["10-01", "12-31"]]

for year = 2000:2022
  for q in quarters
    push!(time_windows,[@sprintf("%d-%s",year,q[1]),@sprintf("%d-%s",year,q[2])])
  end
end

# list of all selected T-wave pairs
floats = DataFrame(x=Float64[], y=Float64[], z0=Float64[], t=DateTime[],
                   Δsp1=Float64[], Δsp2=Float64[], zmax=Float64[], smax=Float64[])

d0 = 360
argoplt = false
rc("font", size=8)
rc("axes", titlesize="medium")
for tw in time_windows
  params = [135,170,15.,45.,0.,2000.,tw[1],tw[2]]
  ntry = 0
  while ntry < 10 
    try
      da = loader.region(params).to_xarray()
      df = DataFrame(x = Array{Float64}(da.LONGITUDE.values),
                     y = Array{Float64}(da.LATITUDE.values),
                     t = Nanosecond.(Array{Int64}(da.TIME.values.astype(np.int64))) + DateTime(1970),
                     p = Array{Float64}(da.PRES.values),
                     #S = Array{Float64}(da.PSAL.values),
                     T = Array{Float64}(da.TEMP.values))

      gdf = groupby(df, [:x,:y,:t])
      
      for (key, subdf) in pairs(gdf)
        ie,je = argmin(abs.(xe.-key.x)),argmin(abs.(ye.-key.y))
        z = gsw.z_from_p.(subdf.p, key.y)
        ik0 = argmin(z)
        ia,ja = argmin(abs.(xa.-key.x)),argmin(abs.(ya.-key.y))
        Taij = Ta[ia,ja,:]
        ka = findlast(.!isnan.(Taij))
        if isnothing(ka)
          continue
        end
        if z[ik0]<-1.9e3 && H[ie,je]>2e3 && -dpa[ka+1]<z[ik0]#sum(isnan.(Ta[ia,ja,:]))==0
          #dists = minimum(SOT.dist.(xks, yks, key.x, key.y))/1e3
          #distn = minimum(SOT.dist.(xkn, ykn, key.x, key.y))/1e3
          #azm = SOT.azimuth(tstation[2],tstation[1],key.y,key.x)-θ0
          if true#at[iamin] < azm < at[iamax] || dists<d0 || distn<d0
            M,C = getMC(dpa[1:ka+1])
            etpT = linear_interpolation(dpa[1:ka+1], M\(C*Taij[1:ka]),extrapolation_bc=Line()) 
            #etpT = linear_interpolation(pa, M\(C*Ta[ia,ja,:]),extrapolation_bc=Line()) 
            ΔT = subdf.T .- etpT.(-z)#etpT.(subdf.p)
            etpw1 = linear_interpolation(ze, wp[ie,je,:,1],extrapolation_bc=Line()) 
            Δsp1 = ΔT.*etpw1.(z)
            kmax = argmax(abs.(Δsp1))
            smax = 1e9*Δsp1[kmax]
            Δsp1 = -1e6*trapz(z,ΔT.*etpw1.(z))
            etpw2 = linear_interpolation(ze, wp[ie,je,:,2],extrapolation_bc=Line()) 
            Δsp2 = -1e6*trapz(z,ΔT.*etpw2.(z))

            #etpn2G1 = linear_interpolation(ze, n2[ie,je,:].*G[ie,je,:,1],extrapolation_bc=Line()) 
            #E = ones(length(z))
            #klow = ze.<z[ik0]
            #Elow = ones(sum(klow))
            #for n = 1:nmode
            #  etpn2Gn = linear_interpolation(ze, n2[ie,je,:].*G[ie,je,:,n],extrapolation_bc=Line()) 
            #  E = [E etpn2Gn.(z)]
            #  Elow = [Elow n2[ie,je,klow].*G[ie,je,klow,n]]
            #end
            #E,Elow=E[:,2:end],Elow[:,2:end]
            #P = inv(iR+E'*iN*E)
            #β = P*E'*iN*ΔT
            #zelow, ΔTlow = ze[klow],Elow*β
            #Δsp1low = 1e6*trapz(zelow,ΔTlow.*etpw1.(zelow))
            #Δsp2low = 1e6*trapz(zelow,ΔTlow.*etpw2.(zelow))
            #e = mean((ΔT .- E*β).^2)^0.5
            if !isnan(Δsp1)
              push!(floats,[key.x,key.y,z[ik0],key.t,Δsp1,Δsp2,z[kmax],smax])
              
              if mod1(size(floats,1),3000)==30 && argoplt
                fig,ax = subplots(1,3,sharey=true)
                ax[1].plot(ΔT,z/1e3,label="data")
                #ax[1].plot((n2[ie,je,:].*G[ie,je,:,1:nmode])*β,ze/1e3,label="model")
                ax[1].legend(frameon=false)
                ax[1].set_xlabel("temperature anomaly (K)")
                ax[1].set_ylabel("z (km)")
                ax[2].plot(1e9*wp[ie,je,:,1],ze/1e3,label = "mode 1")
                ax[2].plot(1e9*wp[ie,je,:,2],ze/1e3,label = "mode 2")
                ax[2].legend(frameon=false)
                ax[2].set_xlabel("mode (s\$~\\mathrm{km}^{-2}~\\mathrm{K}^{-1}\$)")
                ax[3].plot(1e9*ΔT.*etpw1.(z),z/1e3)
                ax[3].plot(1e9*ΔT.*etpw2.(z),z/1e3)
                #ax[3].plot(1e9*ΔTlow.*etpw1.(zelow),zelow/1e3,color="tab:blue",linestyle="--")
                #ax[3].plot(1e9*ΔTlow.*etpw2.(zelow),zelow/1e3,color="tab:orange",linestyle="--")
                ax[3].set_xlabel("projection (s\$~\\mathrm{km}^{-2}\$)")
                fig.suptitle(@sprintf("x=%.1f, y=%.1f, t=%s",key.x,key.y,key.t))
                fig.tight_layout()
                fig.savefig(@sprintf("results/argo/H11_argo_zT%s.pdf",key.t))
              end
            end
          end
        end
      end
      @printf("%s to %s: %d profiles\n",tw[1],tw[2],size(floats,1))
      break
    catch y
      println(y)
    end
  end
  CSV.write(@sprintf("results/argo/japan_%s_%s.csv",time_windows[1][1][1:7],time_windows[end][2][1:7]), floats)
end

include("../src/SOT.jl")
using .SOT, PyCall, Printf, LinearAlgebra, Dates, HDF5, DataFrames, CSV, Statistics, StatsBase
using Random, Distributions, SparseArrays, PyPlot, Interpolations,NCDatasets,HCubature
PyPlot.matplotlib[:rc]("mathtext",fontset="cm")        #computer modern font 
PyPlot.matplotlib[:rc]("font",family="STIXGeneral")
#rc("font", size=10)
#rc("axes", titlesize="large")
rc("font", size=9)
rc("axes", titlesize="medium")

# in situ data
ctdfile = "results/ship/wod/ctds_kuroshio.csv"
argofile = "results/argo/japan_1997-01_2022-12.csv"

# read rg result
h5fl = "results/argo/kuroshio_rg.h5"
xrg = h5read(h5fl, "x")
yrg = h5read(h5fl, "y")
nxrg,nyrg=length(xrg),length(yrg)
nodes = (xrg, yrg)
srcs = [[141.96,38.58],[142.06,38.87]]
tstations = [[166.90986,19.71786],[166.652,19.283]]
evtpos = [142.85,38.10]
nxk = 321
h5fl2 = "results/argo/kuroshio_rg1922.h5"
trg = Millisecond.([h5read(h5fl, "t");h5read(h5fl2, "t")]) .+ DateTime(2000, 1, 1) 
Δsrg = cat(h5read(h5fl, "Δsp1"),h5read(h5fl2, "Δsp1");dims=3)
ntrg = length(trg)

# azimuth and length limits
θ1,θ2,L=-6.5,3.5,3.2

# plot distribution of in situ data
t0,t1,x0,y0,z0 = DateTime(1997,5,6),DateTime(2021,12,31),141,36,-1.9e3
labels = ["shipboard CTDs", "Argo profiles"]
fig,ax=subplots(2,1,figsize=(190/25.4, 190/25.4/2),sharex=true)
c2 = ["k","tab:red"]
for (i,file) in enumerate([ctdfile,argofile])
    floats = CSV.read(file, DataFrame)[:,1:7]
    floats = floats[(floats.z0.<z0),:]
    floats = floats[(floats.y .<=y0) .| (floats.x .>= x0),:]
    floats = floats[(t0 .≤ floats.t .< t1),:]
    floats = floats[.!(isnan.(floats.Δsp1)),:]
    floats = floats[abs.(floats.Δsp1) .< 10,:]
    floats = floats[floats.zmax.>z0,:]
    unique!(floats, [:x,:y,:t])
    sort!(floats, [:t])
    df = combine(groupby(transform(floats, :t => ByRow(year)), :t_year), nrow=>:count)
    edges = DateTime.([df.t_year; max(df.t_year...)+1]) 
    if i==1
        df.count[1] = round(365/240*df.count[1])
        ax[2].set_xlim([edges[1],edges[end]])
    end
    
    label0 = @sprintf("%s",labels[i]) 
    p1 = ax[1].stairs((df.count),edges,color=c2[i],label=label0)
    a0 = SOT.azimuth.(tstations[1][2],tstations[1][1],evtpos[2],evtpos[1])
    floats.a = SOT.azimuth.(tstations[1][2],tstations[1][1],floats.y,floats.x) .- a0
    floats = floats[θ1.<= floats.a .<= θ2,:]
    df = combine(groupby(transform(floats, :t => ByRow(year)), :t_year), nrow=>:count)
    edges = DateTime.([df.t_year; max(df.t_year...)+1]) 
    p2 = ax[2].stairs((df.count),edges,color=c2[i],label=@sprintf("%s",labels[i]))

    ax[i].set_ylim([1,1e4])
    ax[i].set_ylabel("count")
    ax[i].set_yscale("log")
end
ax[1].set_title("\${135^\\circ}\$E to \${170^\\circ}\$E, \${15^\\circ}\$N to \${45^\\circ}\$N")
ax[2].set_title("circular sector to H11")
ax[1].set_title("(a)",loc="left")
ax[2].set_title("(b)",loc="left")
ax[1].legend(frameon=false,loc="upper left")
fig.tight_layout()
fig.savefig("results/argo/point_time.pdf",bbox_inches="tight",dpi=300)

@printf("Now plot quantiles...\n")
p = 0.5:0.5:99.5
α = 0.05
zτh = h5read("results/argo/argotris_H11.h5", "zτ")
zτw = h5read("results/argo/argotris_WAKE.h5", "zτ")
ziP = h5read("results/argo/argotris_ctd.h5", "zs")
println([std(zτh),std(zτw),std(ziP)])
qsc = percentile(ziP,p)
qs4,qs5 = percentile(zτh,p),percentile(zτw,p)
qt = quantile.(Normal(), p/100)
lb,ub = zeros(length(p),3),zeros(length(p),3)
for i = 1:length(p)
    for (j,nt) in enumerate([length(ziP),length(zτh),length(zτw)])
        d = NoncentralT(nt-1,-sqrt(nt)*qt[i])
        lb[i,j],ub[i,j] = -quantile(d, 1-α/2)/sqrt(nt),-quantile(d, α/2)/sqrt(nt)
    end
end
rc("font", size=10)
fig,ax=subplots(3,1,figsize=(4.8,4.8),sharex=true)
c3 = ["#1b9e77","#d95f02","#7570b3"]
ax[1].set_title("H11")
ax[1].set_title("(a)",loc="left")
ax[1].plot(qt,qs4-qt,label="H11",c=c3[1])
ax[1].fill_between(qt, (lb[:,2]-qt), (ub[:,2]-qt), alpha=.2, zorder=3, color=c3[1], linewidth=0)
ax[2].set_title("WAKE")
ax[2].set_title("(b)",loc="left")
ax[2].plot(qt,qs5-qt,label="WAKE",c=c3[2])
ax[2].fill_between(qt, (lb[:,3]-qt), (ub[:,3]-qt), alpha=.2, zorder=3, color=c3[2], linewidth=0)
ax[3].set_title("Shipboard CTD")
ax[3].set_title("(c)",loc="left")
ax[3].plot(qt,qsc-qt,label="shipboard CTD",c=c3[3])
ax[3].fill_between(qt, (lb[:,1]-qt), (ub[:,1]-qt), alpha=.2, zorder=3, color=c3[3], linewidth=0)
#ax.axhline(0,color="black",ls=":",lw=1)
ax[3].set_xlabel("theoretical quantile")
ax[2].set_ylabel("sample quantile \$-\$ theoretical quantile")
#fig.text(0.01, 0.5, "sample quantile \$-\$ theoretical quantile", va="center", rotation="vertical")
#ax.legend(frameon=false)
ax[1].set_ylim([-0.5,0.5])
ax[2].set_ylim([-0.5,0.5])
ax[3].set_xlim([qt[1],qt[end]])
fig.tight_layout()
fig.savefig("results/argo/quantile3_kuroshio.pdf",bbox_inches="tight",dpi=300)

# plot path-mean temperature anomalies
@printf("max %.2e, min %.2e\n",max(Δsrg...),min(Δsrg...))
σs = [0.44,0.44]
itpc0 = SOT.cint(0)
σtrend = 0.037/SOT.meanyear
K = -6
flnms = ["results/argo/point_H11_0d.h5","results/argo/point_WAKE_.6d.h5"]
starts = [DateTime(2007,7,1),DateTime(1996,12,1)]
#tes = [DateTime(2008, 1, 1):Day(10):DateTime(2021, 7, 31),DateTime(1997, 6, 1):Day(10):DateTime(2021, 7, 31)]
fig,ax=subplots(4,1,figsize=(190/25.4, 190/25.4),gridspec_kw=Dict("height_ratios"=>(3,2,3,2)))
for (i,file) in enumerate(flnms)
    xa = h5read(file, "xa")
    xt = h5read(file, "xt")
    xat = h5read(file, "xat")
    Pa = h5read(file, "Pa")
    Pt = h5read(file, "Pt")
    Pat = h5read(file, "Pat")
    ter = h5read(file, "t")
    te = DateTime(2000, 1, 1) .+ Millisecond.(ter)
    ter /= (1000*24*3600)
    tm = (max(ter...)+min(ter...))/2
    ω = 2π/SOT.meanyear
    m = length(ter)
    D = [I(m) Diagonal(ter.-tm) cos.(ω*ter) sin.(ω*ter) cos.(2ω*ter) sin.(2ω*ter)]
    τa,τt,τat = D*xa,D*xt,D*xat
    ea,et,eat = sqrt.(diag(D*Pa*D')),sqrt.(diag(D*Pt*D')),sqrt.(diag(D*Pat*D'))
    e0 = vec(h5read(file, "e0"))
    xk,yk = SOT.findpath(srcs[i], tstations[i], nxk)
    svrt,svra,svrj = 1-mean(diag(Pt)[1:m])/(σs[i]^2*itpc0),1-mean(diag(Pa)[1:m])/(σs[i]^2*itpc0),1-mean(diag(Pat)[1:m])/(σs[i]^2*itpc0)
    @printf("stochastic variance reduction: Twave %.2f, Argo %.2f, joint %.2f\n",svrt,svra,svrj)
    @printf("trend variance reduction: Twave %.2f, Argo %.2f, joint %.2f\n",1-Pt[m+1,m+1]/σtrend^2,1-Pa[m+1,m+1]/σtrend^2,1-Pat[m+1,m+1]/σtrend^2)
    @printf("full variance reduction: Twave %.2f, Argo %.2f, joint %.2f\n",1-(sum(et.^2)/sum(e0.^2)),1-(sum(ea.^2)/sum(e0.^2)),1-(sum(eat.^2)/sum(e0.^2)))
    #println(1-(sum(et.^2)/sum(ea.^2)))
    println(1-sum(diag(Pat)[1:m])/sum(diag(Pa)[1:m]))
    println(1-(diag(Pat)[m+1]/diag(Pa)[m+1]))
    println(1-(sum(eat.^2)/sum(ea.^2)))
    @printf("Argo: trend %.1e K/yr, etrend %.1e K/yr\n",xa[m+1]*SOT.meanyear/K,-2sqrt(diag(Pa)[m+1])*SOT.meanyear/K)
    @printf("Twave: trend %.1e K/yr, etrend %.1e K/yr\n",xt[m+1]*SOT.meanyear/K,-2sqrt(diag(Pt)[m+1])*SOT.meanyear/K)
    @printf("joint: trend %.1e K/yr, etrend %.1e K/yr\n",xat[m+1]*SOT.meanyear/K,-2sqrt(diag(Pat)[m+1])*SOT.meanyear/K)
    
    τrg = zeros(ntrg)
    for j = 1:ntrg
        itprg = interpolate(nodes, Δsrg[:,:,j], Gridded(Linear()))
        τrg[j] = sum(itprg.(xk,yk))*3.2/nxk
    end
    if i == 2
      trgd = Dates.value.(trg .- DateTime(2000,1,1))/(1000*24*3600)
      itprg = interpolate((trgd,), τrg, Gridded(Linear()))
      irg = te .> trg[1]
      ti,tri = te[irg],ter[irg]
      println(ti[argmax(τat[irg].-itprg.(tri))])
    end

    #τa,τt,τat = xa[1:m],xt[1:m],xat[1:m]
    #ea,et,eat = sqrt.(diag(Pa))[1:m],sqrt.(diag(Pt))[1:m],sqrt.(diag(Pat))[1:m]
    va,vt,vat = diag(Pa)[1:m]/K^2,diag(Pt)[1:m]/K^2,diag(Pat)[1:m]/K^2
    ax[2*i-1].set_ylim([-0.3,0.3])
    ax[2*i-1].plot(trg,τrg/K, linewidth=1, linestyle="--",color="tab:blue")
    p1, = ax[2*i-1].plot(te,τa/K, linewidth=1,label="point mode 1 ",color="tab:blue")
    ax[2*i-1].fill_between(te, (τa-2ea)/K, (τa+2ea)/K, alpha=.2, zorder=3, color="tab:blue", linewidth=0)
    p2, = ax[2*i-1].plot(te,τt/K, linewidth=1,label="\$T\$ waves",color="tab:orange")
    ax[2*i-1].fill_between(te, (τt-2et)/K, (τt+2et)/K, alpha=.2, zorder=3, color="tab:orange", linewidth=0)
    p3, = ax[2*i-1].plot(te,τat/K, linewidth=1,label="point and \$\\mathrm{T}\$ waves",color="tab:green")
    ax[2*i-1].fill_between(te, (τat-2eat)/K, (τat+2eat)/K, alpha=.2, zorder=3, color="tab:green", linewidth=0)
    ax[2*i-1].set_ylabel("temperature anomaly (K)")
    ax[2*i-1].set_xlim([te[1],te[end]])

    ax[2*i].set_ylim([0,(σs[i]/K)^2*itpc0])
    #p4 = ax[2*i].axhline(2*σs[i]*sqrt(itpc0)/abs(K),color="tab:red",linewidth=1,zorder=0)
    #p5, = ax[2*i].plot(te,2*e0/abs(K),color="k", linewidth=1, zorder=0)
    ax[2*i].plot(te,va,color="tab:blue", linewidth=1)
    ax[2*i].plot(te,vt,color="tab:orange", linewidth=1)
    ax[2*i].plot(te,vat,color="tab:green", linewidth=1)
    ax[2*i].set_ylabel("variance (\$\\mathrm{K}^2\$)")
    ax[2*i].set_xlim([te[1],te[end]])
    ax[2*i].set_title("($(('a':'z')[2*i]))", loc="left")
    if i==1
      ax[1].set_title("(a)", loc="left")
      ax[1].legend([p1, p2, p3], ["conventional","seismic","combined"], ncol=3, loc="upper center", frameon=false)
    #  ax[2].legend([p4,p5],["stochastic prior","full prior"],ncol=2,loc="upper center",frameon=false)
    end
end
ax[3].set_title("(c)", loc="left")
ax[1].set_title("\$\\theta={316.5^\\circ}\$ to H11")
ax[3].set_title("\$\\theta={317.1^\\circ}\$ to WAKE")
fig.align_ylabels()
fig.tight_layout()
fig.savefig("results/argo/point_path.pdf",bbox_inches="tight",dpi=300)

# plot area-mean temperature anomalies
λx = 100e-3/sqrt(2)
shape = "tri"
if shape=="tri"
    c0(x) = sqrt(2)*x[1]*x[2]*exp(-sqrt(x[1]^2+x[2]^2-2x[1]*x[2]*cosd(x[3]-x[4]))/λx/sqrt(2))*cos(sqrt(x[1]^2+x[2]^2-2x[1]*x[2]*cosd(x[3]-x[4]))/λx/sqrt(2)-π/4)
    c2d = hcubature(c0, [0,0,θ1,θ1], [L,L,θ2,θ2]; rtol=1e-6, atol=1e-6)[1]/(((θ2-θ1)*L/2)^2)
    c0(x) = x[1]*x[2]*exp.(-sqrt(x[1]^2+x[2]^2-2x[1]*x[2]*cosd(x[3]-x[4]))/0.2)
    tc2d = hcubature(c0, [0,0,θ1,θ1], [L,L,θ2,θ2]; rtol=1e-6, atol=1e-6)[1]/(((θ2-θ1)*L/2)^2)
else
    c0(x) = sqrt(2)*exp(-sqrt((x[1]-x[2])^2+(x[3]-x[4])^2)/λx/sqrt(2))*cos(sqrt((x[1]-x[2])^2+(x[3]-x[4])^2)/λx/sqrt(2)-π/4)
    c2d = hcubature(c0, [0,0,π*L*θ1/180,π*L*θ1/180], [L,L,π*L*θ2/180,π*L*θ2/180]; rtol=1e-6, atol=1e-6)[1]/((π*L*(θ2-θ1)/180)^2)
    c0(x) = exp.(-sqrt((x[1]-x[2])^2+(x[3]-x[4])^2)/0.2)
    tc2d = hcubature(c0, [0,0,π*L*θ1/180,π*L*θ1/180], [L,L,π*L*θ2/180,π*L*θ2/180]; rtol=1e-6, atol=1e-6)[1]/((π*L*(θ2-θ1)/180)^2)
end
tc2d *= σtrend^2/SOT.tcint(0)
println(-sqrt(tc2d)*SOT.meanyear/K)
fig,ax=plt.subplots(4,1,figsize=(190/25.4, 190/25.4),gridspec_kw=Dict("height_ratios"=>(3,2,3,2)))
for (i,receiver) in enumerate(["H11","WAKE"])
    println(receiver)
    h5name = @sprintf("results/argo/point_%s_2d%s.h5",receiver,shape)
    xa = h5read(h5name, "xa")
    xt = h5read(h5name, "xt")
    xat = h5read(h5name, "xat")
    Pa = h5read(h5name, "Pa")
    Pt = h5read(h5name, "Pt")
    Pat = h5read(h5name, "Pat")
    ter = h5read(h5name, "t")
    e0 = vec(h5read(h5name, "e0"))
    te = DateTime(2000, 1, 1) .+ Millisecond.(ter)
    ter /= (1000*24*3600)
    m = length(te)
    tm = (max(ter...)+min(ter...))/2
    ω = 2π/SOT.meanyear
    D = [I(m) (ter.-tm) cos.(ω*ter) sin.(ω*ter) cos.(2ω*ter) sin.(2ω*ter)]
    τa,τt,τat = D*xa,D*xt,D*xat
    ea,et,eat = sqrt.(diag(D*Pa*D')),sqrt.(diag(D*Pt*D')),sqrt.(diag(D*Pat*D'))
    svrt,svra,svrj = 1-mean(diag(Pt)[1:m])/(σs[i]^2*c2d),1-mean(diag(Pa)[1:m])/(σs[i]^2*c2d),1-mean(diag(Pat)[1:m])/(σs[i]^2*c2d)
    @printf("stochastic variance reduction: Twave %.2f, Argo %.2f, joint %.2f\n",svrt,svra,svrj)
    @printf("trend variance reduction: Twave %.2f, Argo %.2f, joint %.2f\n",1-Pt[m+1,m+1]/tc2d,1-Pa[m+1,m+1]/tc2d,1-Pat[m+1,m+1]/tc2d)
    @printf("full variance reduction: Twave %.2f, Argo %.2f, joint %.2f\n",1-(sum(et.^2)/sum(e0.^2)),1-(sum(ea.^2)/sum(e0.^2)),1-(sum(eat.^2)/sum(e0.^2)))
    #println(1-(sum(et.^2)/sum(ea.^2)))
    println(1-sum(diag(Pat)[1:m])/sum(diag(Pa)[1:m]))
    println(1-(diag(Pat)[m+1]/diag(Pa)[m+1]))
    println(1-(sum(eat.^2)/sum(ea.^2)))
    @printf("Argo: trend %.2e K/yr, etrend %.1e K/yr\n",xa[m+1]*SOT.meanyear/K,-2sqrt(diag(Pa)[m+1])*SOT.meanyear/K)
    @printf("Twave: trend %.2e K/yr, etrend %.1e K/yr\n",xt[m+1]*SOT.meanyear/K,-2sqrt(diag(Pt)[m+1])*SOT.meanyear/K)
    @printf("joint: trend %.2e K/yr, etrend %.1e K/yr\n",xat[m+1]*SOT.meanyear/K,-2sqrt(diag(Pat)[m+1])*SOT.meanyear/K)
    
    if i==1
      global ter0 = ter
      global τt0 = τt
    else 
      itp = linear_interpolation(ter, τt)
      @printf("H11 vs WAKE std %.2e K\n\n",std(itp.(ter0[1:end-1]).-τt0[1:end-1]))
    end
    
    a0 = SOT.azimuth.(tstations[i][2],tstations[i][1],evtpos[2],evtpos[1])
    arg = SOT.azimuth.(tstations[i][2],tstations[i][1],yrg',xrg) .- a0
    drg = SOT.dist.(tstations[i][1], tstations[i][2], xrg, yrg')/1e6
    aidx = θ1.<= arg .<= θ2#abs.(drg.*sind.(arg)) .<= π*3.2*5/180
    
    τrg = [sum(aidx.*Δsrg[:,:,j])/sum(aidx)*3.2 for j = 1:ntrg]
    
    #τa,τt,τat = xa[1:m],xt[1:m],xat[1:m]
    #ea,et,eat = sqrt.(diag(Pa))[1:m],sqrt.(diag(Pt))[1:m],sqrt.(diag(Pat))[1:m]
    va,vt,vat = diag(Pa)[1:m]/K^2,diag(Pt)[1:m]/K^2,diag(Pat)[1:m]/K^2
    p1, = ax[2*i-1].plot(te,τa/K, linewidth=1,color="tab:blue")
    ax[2*i-1].plot(trg,τrg/K, linewidth=1, linestyle="--",color="tab:blue")
    ax[2*i-1].fill_between(te, (τa-2ea)/K, (τa+2ea)/K, alpha=.2, zorder=3, color="tab:blue", linewidth=0)
    p2, = ax[2*i-1].plot(te,τt/K, linewidth=1,label="\$T\$ waves",color="tab:orange")
    ax[2*i-1].fill_between(te, (τt-2et)/K, (τt+2et)/K, alpha=.2, zorder=3, color="tab:orange", linewidth=0)
    p3, = ax[2*i-1].plot(te,τat/K, linewidth=1,label="point and \$\\mathrm{T}\$ waves",color="tab:green")
    ax[2*i-1].fill_between(te, (τat-2eat)/K, (τat+2eat)/K, alpha=.2, zorder=3, color="tab:green", linewidth=0)
    ax[2*i-1].set_ylabel("temperature anomaly (K)")
    ax[2*i-1].set_xlim([te[1],te[end]])
    ax[2*i-1].set_ylim([-0.3,0.3])
    
    ax[2*i].set_ylim([0,(σs[i]/K)^2*c2d])
    #p4 = ax[2*i].axhline(2*σs[i]*sqrt(c2d)/abs(K),color="tab:red",linewidth=1,zorder=0)
    #p5, = ax[2*i].plot(te,2*e0/abs(K),color="k", linewidth=1, zorder=0)
    ax[2*i].plot(te,va,color="tab:blue", linewidth=1)
    ax[2*i].plot(te,vt,color="tab:orange", linewidth=1)
    ax[2*i].plot(te,vat,color="tab:green", linewidth=1)
    ax[2*i].set_ylabel("variance (\$\\mathrm{K}^2\$)")
    ax[2*i].set_xlim([te[1],te[end]])
    ax[2*i].set_title("($(('a':'z')[2*i]))", loc="left")
    if i==1
      ax[1].set_title("(a)", loc="left")
      ax[1].legend([p1, p2, p3], ["conventional","seismic","combined"], ncol=3, loc="upper center", frameon=false)
      #ax[2].legend([p4,p5],["stochastic prior","full prior"],ncol=2,loc="upper center",frameon=false)
    end
    if shape=="tri"
        ax[2*i-1].set_title(@sprintf("Sector to %s",receiver))
    else
        ax[2*i-1].set_title(@sprintf("Rectangular area to %s",receiver))
    end
end
ax[3].set_title("(c)", loc="left")
fig.align_ylabels()
fig.tight_layout()
fig.savefig(@sprintf("results/argo/point_sbasin%s.pdf",shape),bbox_inches="tight",dpi=300)

# plot pointwise trend and offset
rc("font", size=8.5)
Kp = K/3.2
ywindow,yly = ["2008\$-\$2021","1997\$-\$2021"], ["14-yr","25-yr"]
filename = "data/ssh/dt_global_twosat_phy_l4_20181124_vDT2021.nc"
fl = Dataset(filename, "r")
lon = Array{Float64}(fl["longitude"][:])
lat = Array{Float64}(fl["latitude"][:])
ix,iy = 135 .<= lon .<= 170,15 .<= lat .<= 45
sla = Array{Float64}(replace(fl["sla"][ix,iy,1],(missing=>NaN)))
mask = replace(isnan.(sla),(false=>NaN))
fig,axs=plt.subplots(4,2,figsize=(190/25.4-3, 190/25.4),sharex=true,sharey=true)
xk0, yk0 = SOT.findpath(srcs[1], tstations[1], nxk)
axs[1,1].text(xk0[end]-4,yk0[end]-2,"Wake")
axs[1,1].text(xk0[end]-4,yk0[end]-4,"Island")
axs[1,1].text(xk0[1],yk0[1]+1.5,"Japan")
atrend = 55
for (i,receiver) in enumerate(["H11","WAKE"])
    for j = 1:2
      axs[2*i-1,j].set_aspect(1/cosd(30))
      axs[2*i,j].set_aspect(1/cosd(30))
    end
    h5name = @sprintf("results/argo/kuroshio_mt_%s.h5",receiver)
    gmn,gtd = h5read(h5name, "gmn"),h5read(h5name, "gtd")
    
    im,jm = Tuple(argmin(replace(gtd[:,:,1],(NaN=>0))[xrg .< 150,:]))
    @printf("Strongest trend %.2f pm %.2f mK/yr\n",1e3SOT.meanyear*gtd[im,jm,1]/Kp,-2e3SOT.meanyear*gtd[im,jm,2]/Kp)
    println([minimum(1e3SOT.meanyear*replace(gtd[:,:,2],(NaN=>Inf))/abs(Kp)),maximum(1e3SOT.meanyear*replace(gtd[:,:,2],(NaN=>0))/abs(Kp))])
    println([minimum(replace(gmn[:,:,2],(NaN=>Inf))/abs(Kp)),maximum(replace(gmn[:,:,2],(NaN=>0))/abs(Kp))])
    
    local ax = axs[2*i,1]
    if i==1
      global im1=ax.pcolormesh(xrg,yrg,vmin=-1,vmax=1,gmn[:,:,1]'/Kp,cmap="RdBu_r",shading="auto",rasterized=true)
    else
      ax.pcolormesh(xrg,yrg,vmin=-1,vmax=1,gmn[:,:,1]'/Kp,cmap="RdBu_r",shading="auto",rasterized=true)
    end
    ax.pcolormesh(lon[ix],lat[iy],mask',cmap="Greys",vmin=0,vmax=2, shading="auto",rasterized=true)
    ax.set_title(@sprintf("%s offset",ywindow[i]))
    ax.set_title("($(('a':'z')[4*i-1]))",loc="left")
    ax.set_ylabel("latitude")

    local ax = axs[2*i,2]
    ax.set_title(@sprintf("%s offset uncertainty",yly[i]))
    ax.set_title("($(('a':'z')[4*i]))",loc="left")
    if i==1
      global im2=ax.pcolormesh(xrg,yrg,vmin=0.03,vmax=0.14,gmn[:,:,2]'/abs(Kp),shading="auto",rasterized=true)
    else
      im2=ax.pcolormesh(xrg,yrg,vmin=0.03,vmax=0.14,gmn[:,:,2]'/abs(Kp),shading="auto",rasterized=true)
    end
    ax.pcolormesh(lon[ix],lat[iy],mask',cmap="Greys",vmin=0,vmax=2, shading="auto",rasterized=true)
    
    local ax = axs[2i-1,1]
    if i==1
      global im3=ax.pcolormesh(xrg,yrg,1e3SOT.meanyear*gtd[:,:,1]'/Kp,vmin=-atrend,vmax=atrend,cmap="RdBu_r",shading="auto",rasterized=true)
    else
      ax.pcolormesh(xrg,yrg,1e3SOT.meanyear*gtd[:,:,1]'/Kp,vmin=-atrend,vmax=atrend,cmap="RdBu_r",shading="auto",rasterized=true)
    end
    ax.pcolormesh(lon[ix],lat[iy],mask',cmap="Greys",vmin=0,vmax=2, shading="auto",rasterized=true)
    ax.set_title(@sprintf("%s trend",ywindow[i]))
    ax.set_title("($(('a':'z')[4*i-3]))",loc="left")
    ax.set_ylabel("latitude")
    
    local ax = axs[2i-1,2]
    ax.set_title(@sprintf("%s trend uncertainty",yly[i]))
    ax.set_title("($(('a':'z')[4*i-2]))",loc="left")
    if i==1
      global im4=ax.pcolormesh(xrg,yrg,vmin=4,vmax=15,1e3SOT.meanyear*gtd[:,:,2]'/abs(Kp),shading="auto",rasterized=true)
    else
      im4=ax.pcolormesh(xrg,yrg,vmin=4,vmax=15,1e3SOT.meanyear*gtd[:,:,2]'/abs(Kp),shading="auto",rasterized=true)
    end
    ax.pcolormesh(lon[ix],lat[iy],mask',cmap="Greys",vmin=0,vmax=2, shading="auto",rasterized=true)
    
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
    println([[xp[1] yp[1]];[xp[2] yp[2]]])
    xk1,yk1 = SOT.findpath([xp[2],yp[2]],tstations[i],nxk)
    xk2,yk2 = SOT.findpath([xp[1],yp[1]],tstations[i],nxk)
    axs[2i-1,1].plot(xk1,yk1,c="k",lw=.8)
    axs[2i-1,1].plot(xk2,yk2,c="k",lw=.8)
    axs[2i-1,1].plot(xp[1:2],yp[1:2],c="k",lw=.8)
    
    iv,jv = Tuple(argmax(replace(1 .- (gtd[:,:,6]./gtd[:,:,2]).^2,(NaN=>0))))
    @printf("Strongest variance reduction %.2e at (%.2f,%.2f) \n",1 - (gtd[iv,jv,6]/gtd[iv,jv,2])^2,xrg[iv],yrg[jv])

end
axs[4,1].set_xlabel("longitude")
axs[4,2].set_xlabel("longitude")
fig.tight_layout()
subplots_adjust(wspace=0.,hspace=0.2)
y2,y1,dy = .795,.565,.465
for y = [y1, y1-dy]
  cbaxes = fig.add_axes([0.16, y, 0.08, 0.01]) 
  cbaxes.text(1.4, -2.4,"(K)")
  cbar = fig.colorbar(im1, cax=cbaxes, ticks=[-1,0,1], orientation="horizontal")
  cbaxes = fig.add_axes([0.62, y, 0.08, 0.01]) 
  cbaxes.text(0.18, -2.4,"(K)")
  cbar = fig.colorbar(im2, cax=cbaxes, ticks=[0.03,0.14], orientation="horizontal")
end
for y = [y2, y2-dy]
  cbaxes = fig.add_axes([0.165, y, 0.08, 0.01]) 
  cbaxes.text(atrend+25, -2.4,"(mK/yr)")
  cbar = fig.colorbar(im3, cax=cbaxes, ticks=[-atrend,0,atrend], orientation="horizontal")
  cbaxes = fig.add_axes([0.6, y, 0.08, 0.01]) 
  cbaxes.text(17.5, -2.4,"(mK/yr)")
  cbar = fig.colorbar(im4, cax=cbaxes, ticks=[5,10,15], orientation="horizontal")
end
fig.savefig("results/argo/kuroshio_mt.pdf",bbox_inches="tight",dpi=300)

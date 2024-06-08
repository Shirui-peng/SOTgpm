include("../src/SOT.jl")
using .SOT, Printf, LinearAlgebra, Dates, HDF5, DataFrames, CSV, Statistics, StatsBase

# Argo data file location
argofile = "results/argo/argo_H11.csv"

# read Argo profiles 
floats = CSV.read(argofile, DataFrame)

# path-path trend auto-covariance 
itc0 = SOT.tcint(0)

# path-path trend prior (mK/yr)
σtrend = 0.037

# annual cycle prior (s)
σannual = 0.1

# initial guess for MLE
x0 = log.([60,100/sqrt(2),0.5,0.5])

# define negative log-likelihood
f(x) = -1 .* SOT.loglikelihoodargo(x, floats.Δsp1, floats.t, floats.x, floats.y;
                                   σtrend,σannual,itc0,full=true)

# minimize negative log-likelihood 
xk,Hk = SOT.BFGS(f,x0,1e-3,100; xrtol=1e-6)

# print results
@printf("MLE: %s\n",exp.(xk))
@printf("1-σ lb: %s\n ub: %s\n",exp.(xk.-sqrt.(diag(Hk))),exp.(xk.+sqrt.(diag(Hk))))
@printf("2-σ lb: %s\n ub: %s\n",exp.(xk.-2*sqrt.(diag(Hk))),exp.(xk.+2*sqrt.(diag(Hk))))

"""
Example output:

H11
αk = 1.0
k = 13, ek = 3.52e+00
xk = [3.9791655776081702, 4.2511023550552, -0.6414417745434671, -0.8243907478787327, -4.9071303379307025]
.........
αk = 1.0
k = 25, ek = 5.65e-03
xk = [3.9797318474363332, 4.251497535944772, -0.6410938234268898, -0.8244911684036159, -8.962930578665425]

WAKE
αk = 1.0
k = 20, ek = 6.80e-01
xk = [4.091783531544084, 4.287075916602848, -0.5979968956021895, -0.8264226048225342, -5.696133584449357]
...
αk = 1.0
k = 23, ek = 6.99e-02
xk = [4.091659815871754, 4.287090918577277, -0.5980263206387401, -0.8262554556052942, -6.700397665054765]
"""
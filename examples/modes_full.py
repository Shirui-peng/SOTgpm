#!export OMP_NUM_THREADS=1
import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as de
import pandas as pd
from scipy.interpolate import interp1d
import xarray as xr
import time

# 1d extrapolation
def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return np.array(list(map(pointwise, np.array(xs))))

    return ufunclike

def cdif(x):
    return np.insert(x[2:]-x[:-2],[0,-1],[x[1]-x[0],x[-1]-x[-2]])
def grad(x,y):
    return cdif(y)/cdif(x)

lat = None
earthradius = 6371e3 
model = 'ofes'
if lat==None:
    da = xr.open_dataarray(f'data/{model}/sig1xz_Nias_H08_2005.nc')
else:
    if model == 'ofes':
        da = xr.open_dataarray(f'data/{model}/sig1xz_Nias_H08_2005_{lat}S.nc')
    else:
        ds2 = xr.open_dataset(f'data/{model}/sigma1_2005_tile2_mean_wof.nc')
        ds5 = xr.open_dataset(f'data/{model}/sigma1_2005_mean_wof.nc')
        da = xr.concat([ds2,ds5],dim="lon").sigma1.interp(lat=[lat])
        da = da.assign_coords({"x": (da.lon - 40)*earthradius*np.cos(lat*np.pi/180)*np.pi/180})
        lidx = np.where((da.lon>=35)&(da.lon<=105))[0]
        da = da.isel(lon=lidx,lat=0)

n2 = np.full(da.values.shape,np.nan)
if lat==None:
    xidx = np.where((da.x>=0)&(da.x<=2900e3))[0]
else:
    if model == 'ofes':
        xidx = np.where((da.x>=0)&(da.x<=6300e3))[0]
    else:
        xidx = np.where((da.lon>=35)&(da.lon<=105))[0]

for i in xidx:
    n2[i,:] = -9.8*grad(da.z.values,da.values[i,:])/1e3

ds = xr.Dataset(
    data_vars=dict(
        n2 = (["x", "z"], n2),),
    coords=dict(
        x = da.x.values,
        z = da.z.values,),
        attrs=dict(description=f"mean stratification from {model} 2005 daily"),
)

#ds.to_netcdf(f"data/{model}/n2_Nias_H08_2005.nc")

# grid points
n = 256

# gravity acceleration
ga = 9.8

dtype = np.complex128

nmode = 15
zp = ds.z.values
xk = ds.x.values
nxk,nzk = len(xk),len(zp)
hnk = np.full((nmode,nxk,nzk),np.nan)
cnk = np.full((nmode,nxk),np.nan)
Hk = np.full((nxk),np.nan)
#xidx = np.argwhere((xk>=0)&(xk<=np.max(xk)))

# set up basis
zcoord = de.Coordinate('z')
dist = de.Distributor(zcoord, dtype=dtype)

t1 = time.time()
print(zp[-3:])
mean = True
for i in xidx:
    print(f'x = {xk[i]/1e3} km, t1 = {(time.time()-t1)/60:.2f} min',flush=True)
    if mean:
        n2p = np.nanmean(ds.n2.values,axis=0)
    else:
        n2p = ds.n2.values[i].flatten()
    k = np.argwhere(np.isnan(n2p)).max()+1
    if k+1<nzk:
        H = -zp[k]
        Hk[i] = H
        print(H)
        zbasis = de.Chebyshev(zcoord, n, bounds=(-H, 0))
        
        # generate sound slowness field
        h = dist.Field(name='h', bases=zbasis)
        tau_1 = dist.Field(name='tau_1')
        tau_2 = dist.Field(name='tau_2')
        c = dist.Field(name='c')
        n2 = dist.Field(name='n2', bases=zbasis)
        z = dist.local_grid(zbasis)
        
        # interpolated sound speed profile from data
        knots = zp
        itn2 = interp1d(zp, n2p)
        etn2 = extrap1d(itn2)
        n2['g'] = etn2(z)
        
        # Substitutions
        dz = lambda A: de.Differentiate(A, zcoord)
        lift_basis = zbasis.derivative_basis(1)
        lift = lambda A: de.Lift(A, lift_basis, -1)
        hz = dz(h) + lift(tau_1) # First-order reduction
        
        # Problem
        problem = de.EVP([h, tau_1, tau_2], eigenvalue=c, namespace=locals())
        problem.add_equation("c*n2*h + dz(hz) + lift(tau_2) = 0")
        problem.add_equation("h(z=-H) = 0")
        problem.add_equation("h(z=0) = 0")
        #problem.add_equation("c*ga*h(z=0) - dz(h)(z=0)  = 0")
        
        # Solve
        solver = problem.build_solver()
        solver.solve_dense(solver.subproblems[0])
        evals = np.sort(solver.eigenvalues)
        ievls = np.argsort(solver.eigenvalues)[evals>0][:]
        evals = evals[evals>0][:]
        
        knots = z
        for m, idx in enumerate(ievls[:nmode]):
            solver.set_state(idx, solver.subsystems[0])
            hn = ((h / np.sqrt(de.Integrate(n2*h**2, zcoord)/H)).evaluate()['g']).real
            #hn = hn*np.sign(hn[np.argmax(np.abs(hn))])
            ithn = interp1d(z, hn)
            ethn = extrap1d(ithn)
            hnk[m,i] = ethn(zp)
            hnk[m,i] = hnk[m,i]*np.sign(hnk[m,i,-2])
            hnk[m,i,zp<-H] = np.nan
            cnk[m,i] = 1/np.abs(evals[m])**.5
    if mean:
        dst = xr.Dataset(
            data_vars=dict(
                cn=(["n"], cnk[:,i]),
                hn=(["n","z"], hnk[:,i]), ),
            coords=dict(
                n = 1+np.arange(nmode),
                z = zp,),
            attrs=dict(description=f"first {nmode} dynamical modes"),)
        if lat == None:
            dst.to_netcdf(f"results/{model}/modes_mean_{model}.nc")
        else:
            dst.to_netcdf(f"results/{model}/modes_mean_{model}_{lat:.0f}lat.nc")
        break

if not mean:
    dst = xr.Dataset(
        data_vars=dict(
            H = (["x"], Hk),
            cn=(["n","x"], cnk),
            hn=(["n","x", "z"], hnk), ),
        coords=dict(
            n = 1+np.arange(nmode),
            x = xk,
            z = zp,),
        attrs=dict(description=f"first {nmode} dynamical modes"),)
    if lat == None:
        dst.to_netcdf(f"results/{model}/modes_full_{model}.nc")
    else:
        dst.to_netcdf(f"results/{model}/modes_full_{model}_{lat:.0f}lat.nc")

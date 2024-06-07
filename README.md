# SOTgpm.jl
The code in this package processes seismic, acoustic, and in situ data to infer temperature change in the ocean. The principal functions for the seismic and acoustic data are variants of those in [SOT](https://github.com/joernc/SOT):
1. Obtain an earthquake catalog from ISC.
2. Enhance the catalog through template matching at reference stations.
3. Find P-wave pairs at the reference stations.
4. Cross-correlate the waveforms of the identified P-wave pairs at the T-wave station. 
5. Invert the measured travel time changes for travel time anomalies and apply a cycle skipping correction.
With Gaussian process modeling, parameters in prior covariances are determined by a maximum likelihood estimator. Other functions includes joint inversion with in situ data.

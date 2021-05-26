## Jupyter notebook for fitting Confinement and Shear Neutron Reflectivity Data

The following material describes the software used for Neutron Reflectivity data fitting, for the particular experiments descibed on the paper (to be cited).

## Schematic representation of the sample

## Using RefNX software for modelling data

Refnx is a powerful Python package  that has been used to model and fit the Neutron Reflectometry data obtained with the Confinement Cell Nellie [[1]](#1) and the Confinement and Shear cell Ofelia [[2]](#2).

<p align="center">
<img src="https://user-images.githubusercontent.com/55095328/119741835-cc744880-be86-11eb-9b7c-554517f9eafb.png" width = "800">
</p>

In our experiments, a Mixed Reflectivity model has been used following the provided guidelines for fitting the data using a Jupyter Notebook script from RefNX. Additional information on how to operate with RefNX can be found elsewhere [[3]](#3).  Here, a brief description of a fitted reflectivity experiment performed using Ofelia/Nellie is presented:

1. Packages needed to load the fitting procedure

```python
%matplotlib inline
import sys
import matplotlib.pyplot as plt
import numpy as np
import os.path
import refnx, scipy
# the analysis module contains the curvefitting engine
from refnx.analysis import CurveFitter, Objective, Parameter, GlobalObjective, Transform
# the reflect module contains functionality relevant to reflectometry
from refnx.reflect import SLD, ReflectModel, Structure, MixedReflectModel
# the ReflectDataset object will contain the data
from refnx.dataset import ReflectDataset
# version numbers used in this analysis
refnx.version.version, scipy.version.version
```

2. Data paths

```python
#Full path, where the data is located
pth = '/data'

#Data files: in our case, 4 columns data have been used, 
#q, R (Reflectivity), error_R, error_q
#Silicon block
data_d2oblock = ReflectDataset(os.path.join(pth,'D2O_block.mft'))
data_d2oblock.name = "d2osiblock"

#Sample in D2O contrast
data_d2osample = ReflectDataset(os.path.join(pth,'Sample_in_D2O.mft'))
data_d2osample.name = "d2oSilanessample"

#Sample in H2O contrast
data_h2osample = ReflectDataset(os.path.join(pth,'Sample_in_H2O.mft'))
data_h2osample.name = "h2oSilanessample"

#Polysterene characterization
data_hPS = ReflectDataset(os.path.join(pth,'Polysterene_Melinex.mft'))
data_hPS.name = "hPS"

#Depending on the experiments, this list will be longer: several pressures, several shear speeds
#Pressure 1 bar, sample confined, static
data_d2osampleconfined = ReflectDataset(os.path.join(pth,'P1_sample_confined.mft'))
data_d2osampleconfined.name = "P1_sample_confined"

#Pressure 2 bar, sample confined, static
#data_P2_sample_confined = ReflectDataset(os.path.join(pth,'P2_sample_confined.mft'))
#data_P2_sample_confined.name = "P2_sample_confined"

#Pressure 1 bar, sample confined, v = 100 nm/s
#data_P1_sample_confined = ReflectDataset(os.path.join(pth,'P1_sample_confined.mft'))
#data_P1_sample_confined_v1.name = "P1_sample_confined_v1"

#Pressure 2 bar, sample confined, v = 100 nm/s
#data_P2_sample_confined = ReflectDataset(os.path.join(pth,'P2_sample_confined.mft'))
#data_P2_sample_confined_v1.name = "P2_sample_confined_v1"
```

3. SLD definitions

```python

#SLD of Silicon block layers, D2O, H2O and air (theoretical values)
si = SLD(2.07 + 0j)
sio2 = SLD(3.47 + 0j)
d2o = SLD(6.36 + 0j)
h2o = SLD(-0.56 + 0j)
air = SLD(0.0 + 0j)

#SLD of Silanes (theoretical value)
silanes = SLD(0.7+0J)
silanes.real.setp(vary=False, bounds=(0.6, 0.8))
silanes.real.name='silanes SLD'

#SLD of the unconfined sample in D2O. Use "vary = False" if it is unknown.
sldsample_d2o = Parameter(5.55583, 'sld sample d2o')
sldsample_d2o.setp(vary = False, bounds = (5.2, 6.2))

#SLD of the unconfined sample in H2O. Use "vary = False" if it is unknown.
sldsample_h2o = Parameter(5.82485, 'sld sample h2o')
sldsample_h2o.setp(vary = False, bounds = (3.7, 5.9))

#SLD of the confined sample in D2O (in this case, the same as D2O unconfined).
#sample_nopocket refers to the sample properly confined
#sample_pocket refers to the sample that has been confined with dust or any other contamination
#For this experiment, both are considered equal, as not material variation takes place
Sample_nopocket = SLD(sldsample_d2o)
Sample_pocket = SLD(sldsample_d2o)

#SLD for the backing layer of the pockets contamination
pockets = SLD(4.94199 + 0J)
pockets.real.setp(vary=False, bounds=(4.0, 5.0))
pockets.real.name='pocket SLD 1'

#SLD for the backing layer without pockets contamination
nopockets = SLD(2.61373 + 0J)
nopockets.real.setp(vary=False, bounds=(2.2, 3.0))
nopockets.real.name='no pocket SLD 1'

#SLD of the (hydrogenated) polysterene (theoretical value)
hps = SLD(1.412 +0J)
hps.real.setp(vary=False, bounds= (1.2, 2.8))
hps.real.name='sld hps'

#SLD Melinex (already obtained)
Melinex = SLD(2.53323 + 0J)
Melinex.real.setp(vary=False, bounds=(2.4, 2.8))
Melinex.real.name='silanes SLD'
```

4. Modelling the sample using layers ("snabs") and Mixed Reflectivity. Depending on the data available, several fittings can be made. For example, if data for the unconfined sample, substrate or polysterene is provided, the Mixed Reflectivity will be easier using the results from those fits.
```python
#Si block D2O

sio2_slab = sio2(13.8936, 2.19514)
sio2_slab.thick.setp(vary=False, bounds=(13.8, 14))
sio2_slab.thick.name = 'sio2 thickness'
sio2_slab.rough.setp(vary=False, bounds=(1.9, 8))
sio2_slab.rough.name = name='sio2 roughness'
sio2_slab.vfsolv.setp(0.0, vary=False, bounds=(0.001, 0.4))
sio2_slab.vfsolv.name = 'sio2 solvation'

#Silanes (if the block is silanized)
silanes_slab = silanes(23.4608 ,10.8417)
silanes_slab.thick.setp(vary=False, bounds=(22, 24))
silanes_slab.thick.name = 'silanes thickness'
silanes_slab.rough.setp(vary=False, bounds=(1, 12))
silanes_slab.rough.name = name='silanes/sio2 roughness'
silanes_slab.vfsolv.setp(0.188672, vary=False, bounds=(0.01, 0.3))
silanes_slab.vfsolv.name = 'silanes solvation'
```
<!--- 
solv_roughness1 = Parameter(1.762616, 'solvent roughness d2o')
solv_roughness1.setp(vary=False, bounds=(1, 10))

#Structure to be fitted
s_d2oblock = si | sio2_slab | silanes_slab | d2o(0, solv_roughness1)
s_d2oblock.contract = 1.5

model_d2oblock = ReflectModel(s_d2oblock,scale=0.979909,dq_type = 'pointwise')
model_d2oblock.scale.setp(bounds=(0.9, 1.4), vary=False)
model_d2oblock.bkg.setp(8.17618e-07,bounds=(4e-8, 6e-6), vary=False)

objective_d2oblock = Objective(model_d2oblock, data_d2oblock,use_weights = False)
-->

```python
#Si block with sample and D2O

sampled2o_slab = sampled2o(441.043, 6.17033)
sampled2o_slab.thick.setp(vary=False, bounds=(400, 500))
sampled2o_slab.thick.name = 'sample thickness d2o'
sampled2o_slab.rough.setp(vary=False, bounds=(0.001, 25))
sampled2o_slab.rough.name = name='sample/silanes roughness d2o'
sampled2o_slab.vfsolv.setp(0.924083, vary=False, bounds=(0.6, 0.95))
sampled2o_slab.vfsolv.name = 'sample solvation d2o'
```
<!--- 
solv_roughness2 = Parameter(194.54, 'sample/solvent roughness d2o')
solv_roughness2.setp(vary=False, bounds=(100, 300))

s_d2osample = si | sio2_slab | silanes_slab | sampled2o_slab | d2o(0, solv_roughness2)
s_d2osample.contract = 1.5

model_d2osample = ReflectModel(s_d2osample,scale=1.0655,dq_type = 'pointwise')
model_d2osample.scale.setp(bounds=(0.9, 1.4), vary=False)
model_d2osample.bkg.setp(2.32067e-06,bounds=(4e-12, 1e-5), vary=False)
#model_d2osample.dq.setp(4.6689,bounds=(1, 12), vary=False)

objective_d2osample = Objective(model_d2osample, data_d2osample,use_weights = False)
-->

```python
#Si block with sample and H2O

sampleh2o_slab = sampleh2o(441.043, 4.07343)
sampleh2o_slab.thick.setp(vary=False, bounds=(10, 900))
sampleh2o_slab.thick.constraint = sampled2o_slab.thick
sampleh2o_slab.thick.name = 'sample thickness d2o'
sampleh2o_slab.rough.setp(vary=False, bounds=(0.001, 30))
sampleh2o_slab.rough.constraint = sampled2o_slab.rough
sampleh2o_slab.vfsolv.setp(0.924083, vary=False, bounds=(0.01, 1.0))
sampleh2o_slab.vfsolv.constraint = sampled2o_slab.vfsolv
```
<!--- 
#solv_roughness2 = Parameter(97.4674 , 'sample/solvent roughness')
#solv_roughness2.setp(vary=False, bounds=(0.01, 200))

s_h2osample = si | sio2_slab | silanes_slab | sampleh2o_slab |h2o(0, solv_roughness2)

model_h2osample = ReflectModel(s_h2osample,scale=1.0655, bkg=2.03462e-06, dq_type ='pointwise')
model_h2osample.scale.setp(bounds=(0.9, 1.4), vary=False)
model_h2osample.bkg.setp(bounds=(4e-9, 6e-6), vary=False)
#model_h2osample.dq.setp(bounds=(4, 6), vary=False)

objective_h2osample = Objective(model_h2osample, data_h2osample,use_weights = False)

global_objective = GlobalObjective([objective_d2oblock,objective_d2osample,objective_h2osample])
-->

```python
#hPS

hPS_slab1 = hps(32.9551, 6.44189) 
hPS_slab1.thick.setp(vary=False, bounds=(20, 50))
hPS_slab1.thick.name = 'hPS thickness'
hPS_slab1.rough.setp(vary=False, bounds=(4, 8))
hPS_slab1.rough.name = name='hPS roughness'
hPS_slab1.vfsolv.setp(0.0, vary=False, bounds=(0.0001, 0.4))
hPS_slab1.vfsolv.name = 'hPS solvation'
```
<!--- 
solv_roughness3 = Parameter(22.7224, 'hPS/Melinex roughness')
solv_roughness3.setp(vary=False, bounds=(1, 29))

s_hPS = air | hPS_slab1 | Melinex (0, solv_roughness3)

model_hPS = ReflectModel(s_hPS,scale=1.01611, bkg=1.18178e-06, dq_type ='pointwise')
model_hPS.scale.setp(bounds=(0.9, 1.4), vary=False)
model_hPS.bkg.setp(bounds=(4e-9, 6e-6), vary=False)
#model_hPS.dq.setp(bounds=(4, 6), vary=False)

objective_hPS = Objective(model_hPS, data_hPS,use_weights = False)
-->

```python
#Confined sample 1 Bar
#Confined experiments can be better fitted with a Mixed Reflectivity model. For this model, 
#the layer where the variation happens is the sample layer. One is defined without pockets 
#(the pure confined sample), the other with pockets (the sample with D2O and dust confined 
#between the substrate and the plastic membrane).

samplenopocket_slab = samplenopocket(23.9481, 11.7125)
samplenopocket_slab.thick.setp(vary=False, bounds=(15, 155))
samplenopocket_slab.thick.name = 'sample no pocket thickness 1bar'
samplenopocket_slab.rough.setp(vary=False, bounds=(1, 30))
samplenopocket_slab.rough.name = name='sample/silanes no pocket r 1bar'
samplenopocket_slab.vfsolv.setp(0.306441, vary=False, bounds=(0.01, 1.0))
samplenopocket_slab.vfsolv.name = 'sample no pocket solvation'

samplepocket_slab = samplepocket(93.4844 , 10.9347)
samplepocket_slab.thick.setp(vary=False, bounds=(1, 140))
samplepocket_slab.thick.name = 'sample pocket thickness'
samplepocket_slab.rough.setp(vary=False, bounds=(0.001, 30))
samplepocket_slab.rough.constraint = samplenopocket_slab.rough
samplepocket_slab.vfsolv.setp(0.325543, vary=False, bounds=(0.01, 1.0))
samplepocket_slab.vfsolv.name = 'sample pocket solvation'

solv_roughness1 = Parameter(12.3995, 'No pocket roughness')
solv_roughness1.setp(vary=False, bounds=(3, 20))

solv_roughness2 = Parameter(22.6623, 'Pocket roughness')
solv_roughness2.setp(vary=False, bounds=(17, 30))

#The Mixed Reflectivity model needs a scale value for each reflectivity structure. The scale
#can be used as the "proportion" of each reflectivity, when it is properly normalized.
scale1 = Parameter(0.6039300181837898, 'Scale No pockets')
scale1.setp(vary=False, bounds=(0.6,0.72))
scale2 = Parameter(0.007841735744284222, 'Scale pockets')
scale2.setp(vary=False, bounds=(0.001, 1))

#These are the two structures used for the Mixed Reflectivity
s_d2osampleconfined1 = si | sio2_slab | silanes_slab | samplenopocket_slab | hPS_slab | nopockets(0, solv_roughness1)
s_d2osampleconfined2 = si | sio2_slab | silanes_slab | samplepocket_slab | pockets(0, solv_roughness2)
s_d2osampleconfined1.solvent = SLD(6.36 + 0j)
s_d2osampleconfined2.solvent = SLD(6.36 + 0j)

s_d2osampleconfined = [s_d2osampleconfined1,s_d2osampleconfined2]
scale = [scale1,scale2]

#Preparation of the model to be fitted
model_d2osampleconfined = MixedReflectModel(s_d2osampleconfined, scale, dq_type ='pointwise')
model_d2osampleconfined.bkg.setp(4.5921055645868026e-07,bounds=(1e-7, 9e-5), vary=False)

#Objective function to be fitted
objective_d2osampleconfined = Objective(model_d2osampleconfined, data_d2osampleconfined,use_weights = False)

```
5. Fitting procedure and results
```python
fitter = CurveFitter(objective_d2osampleconfined)
fitter.fit('differential_evolution');
print(objective_d2osampleconfined)
```
6. Refining through MCMC: once the values are obtained from the previous step, a better estimation of the errors (the behaviour of the parameters around the obtained value in comparison with each parameter) on the fitted parameters can be obtained by making a second fit, this time, using the Affine Invariant Markov chain Monte Carlo (MCMC) Ensemble sampler for Bayesian parameter estimation.

```python
fitter = CurveFitter(objective_d2osampleconfined, nwalkers=200)
np.random.seed(6)
fitter.initialise('jitter')
fitter.reset()
fitter.sample(400, random_state=1);
```


7. Results from MCMC 
```python
process_chain(objective_d2osampleconfined, fitter.chain);
print(objective_d2osampleconfined)
objective_d2osampleconfined.corner();
```

8. Plot fit and data
```python 
# Data
plt.rcParams['figure.figsize'] = [15, 12]
plt.errorbar(data_d2osampleconfined.x, data_d2osampleconfined.y, data_d2osampleconfined.y_err,data_d2osampleconfined.x_err,
label='$\mathregular{P=1\ bar}$', ms=4, marker='o',lw=0, elinewidth=1, color='b')

# Fit
plt.plot(data_d2osampleconfined.x, objective_d2osampleconfined.generative(), color='r', zorder=20, linewidth=4)

ax = plt.gca()
ax.set_xscale("log", nonpositive='clip')
ax.set_yscale("log", nonpositive='clip')

plt.legend()
plt.yscale('log')
plt.ylabel('Reflectivity')
plt.xlabel('Q /$\AA^{-1}$')
plt.ylim(0.5*1e-6, 2);
plt.xlim(0.004, 0.2)



```
9. SLD profile
```python 
#SLD

plt.plot(*s_d2omucinsconfined1.sld_profile(),label='$\mathregular{P=1\ bar}$')
plt.plot(*s_d2omucinsconfined2.sld_profile(),label='$\mathregular{P=2\ bar}$')
plt.ylim(-1, 7);
plt.xlim(-40, 800)
plt.legend()
plt.ylabel('SLD /$10^{-6} \AA^{-2}$')
plt.xlabel('distance / $\AA$');

```
## References

<a id="1">[1]</a>  Nelson, A. R., & Prescott, S. W. (2019). refnx: neutron and X-ray reflectometry analysis in Python. Journal of applied crystallography, 52(1), 193-200.

<a id="2">[2]</a>
Nordforsk Project: Neutron scattering of confined and sheared thin soft films. Retrieved from https://www.nordforsk.org/projects/neutron-scattering-confined-and-sheared-thin-soft-films (2021)

<a id="3">[3]</a>  RefNX software.  Retrieved from https://refnx.readthedocs.io/en/latest/#

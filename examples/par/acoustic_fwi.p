# time parameters, default are usually good
courant=0.6
dt=-1
fmax=30
resampling=sinc
sinc_half_length=11
sub=-1

# sources and receivers geometry
ns=10
nr=401
sx0=2
sz0=0.1
rx0=2
rz0=0.1
rxinc=0.01
rzinc=0
sxinc=0.4
szinc=0

# receiver type
seismotype=0

# boundary parameters
bc_top=1
bc_bottom=2
bc_left=2
bc_right=2
taper_top=0
taper_left=50
taper_right=50
taper_bottom=50
taper_strength=0.08

# inversion parameters
niter=10
nlsolver=lbfgs
max_trial=5

# model bounds
vpmin=1.0
vpmax=3.0
rhomin=1.0 
rhomax=2.5

# miscallenous
device=0
nthreads=24
verbose=1
format=0
# time parameters
courant=0.6
dt=-1
fmax=30
resampling=sinc
sinc_half_length=11
sub=20

# sources and receivers geometry
ns=1
nr=101
sx0=1
sz0=0.1
rx0=0.5
rz0=0.2
rxinc=0.01
rzinc=0


# source mechanism and receiver type
seismotype=0

# boundary parameters
bc_top=1
bc_bottom=2
bc_left=2
bc_right=2
taper_top=0
taper_left=30
taper_right=30
taper_bottom=30
taper_strength=0.05

# model bounds
vpmin=0.2
vpmax=8
rhomin=0 
rhomax=8

# miscallenous
device=0
nthreads=24
verbose=3
format=0
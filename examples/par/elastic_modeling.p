# time parameters
courant=0.6
dt=-1
fmax=10
resampling=sinc
sinc_half_length=11
sub=0

# sources and receivers geometry
srcoord=./par/srcoord.txt

# source mechanism and receiver type
mt=1
fangle=0
mxx=0
mzz=0
mxz=1
seismotype=1
gl=0

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
vsmin=0.1
vsmax=5
rhomin=0 
rhomax=8
deltamin=-0.5
deltamax=1
epsilonmin=-0.5
epsilonmax=1

# miscallenous
device=0
nthreads=24
verbose=3
format=0
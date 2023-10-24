##### The following parameters must always be provided in command line and not in a parameter file
# source: single source time function with a single time axis or 2D/3D time function with time, source index, and component (for vector or moment tensor sources) axes
# model: 3D model with depth, X, and parameter axes
# data: 2D/3D data with Time, receiver-source concatenated, and component axes
# output
# wavefield: 3D/4D with Depth, X, component (if many), and Time axes 
# ioutput: prefix for the intermediate outputs of FWI

#############################################
######## time propagation parameters ########
#############################################

# time sampling (in seconds) of the wave propagation. Possible values:
# 0: sampling will be calculated based on Courant number in the least stiff manner.
# -1: sampling will be calculated based on Courant number such as the propagation sampling is an integer multiple of the original source time function sampling.
# val > 0: sampling will be exactly 'val'. It is the responsibility of the user to ensure time stability.
dt=-1

# Courant number for the CFL condition. It may be re-calculated internally based on the choice of 'dt'.
courant=0.6

# maximum frequency (in Hertz) to be propagated. This parameter has no effect on the propagation. It is used for information only to calculate the number of grid points per shortest wavelength.
fmax=10.0

# Type of resampling in time of the source time function and data between original and propagation sampling rates. Possible options 'sinc' and 'linear'.
resampling=sinc

# half-length of the filter (in number of points) used for 'sinc' interpolation
sinc_half_length=11

# subsampling of the wavefield (integer number) when the wavefield is saved or used to compute FWI gradient. Possible values:
# 0: wavefield won't be saved
# val > 0: wavefield (and FWI gradient) will be saved (gradient computed) every 'val' propagation time steps.
# -1: the subsampling 'val' will be calculated internally to remain close but inferior to the Nyquist sampling of the original source time function.
sub=0


#############################################
####### sources-receivers geometry ##########
#############################################

# for arbitrary geometry, provide a file describing source-receiver locations (see example in 'srcoord.txt').
srcoord=srcoord.txt

# Alternatively, a regular line geometry can be provided
# Number of sources
ns=10

# Number of receivers for each source
nr=10

# coordinates of the first source and first receiver (in km)
sx0=0.0
sz0=0.1
rx0=0.0
rz0=0.2

# increment in x and z (in km)
sxinc=0.1
szinc=0.0
rxinc=0.01
rzinc=0.0

# first DAS channel dip (in radians from x-axis clockwise) and increment
rdip0=1.5708
rdipinc=0.0

#############################################
### source mechanism and receiver type ######
#############################################

# source mechanism parameters are relevant for elastic media

# Moment tensor source (1) or vector force (0)
mt=1

# Moment tensor time-independent components
mxx=0
mzz=0
mxz=1

# dipping angle of the vector force (in radians from x-axis clockwise). '0' means horizontal force. 
fangle=0.0

# gauge length for DAS modeling (in km). Possible values:
# 0: no DAS modeling
# val > 0: the receivers will be assumed to be DAS channels with gauge length = val.
gl=0.008

# receiver type (integer). Possible values: 0 or 1 for particle displacement/velocity for geophones, or strain/strain rate for DAS, or hydraulic pressure/time derivative of pressure for hydrophones.
seismotype=1

#############################################
######### boundary conditions ###############
#############################################

# boundary condition. Possible values: 0 (do nothing), 1 (free-surface), 2 (locally absorbing)
bc_top=1
bc_bottom=2
bc_left=2
bc_right=2

# boundary taper (in number of grid points) to be typically used in conjunction with locally absorbing boundary to improve the absorption efficacy
# the model will not be extended spatially, so it is the responsibility of the user to provide a large enough model. 
taper_top=0
taper_left=30
taper_right=30
taper_bottom=30

# efficacy of the cosine-squared tapering zone. '0' means no tapering at all. Larger values attenuate short wavelengths but reflect long wavelengths.
# Typical values are between 0.05 and 0.1
taper_strength=0.05

# free-surface stiffness for acoustic media. Value > 1.
free_surface_stiffness=1.05

#############################################
############ model bounds ###################
#############################################

# model bounds are used for both modeling and inversion
# parameters are used only when relevant
# velocity bounds are in (km/s), density bounds in (g/cc)
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

#############################################
######### inversion parameters ##############
#############################################

# maximum number of iterations
niter=10

# maximum number of trials per iteration
max_trial=5

# convergence-rate threshold to stop the inversion
threshold=0.0

# solver to use: 'lbfgs', 'bfgs', 'nlcg', 'nlsd'
nlsolver=lbfgs

# line search to use: 'weak_wolfe', 'regular_wolfe', 'strong_wolfe'
lsearch=regular_wolfe

# mask to be applied to the gradient. Must have the same geometry as the model. Provide the file name.
mask=none

# weights to be applied to the data. Must have the same size as the data. Provide the file name.
weights=none

# time filter to be applied to the data (as part of the inversion). Provide the file name.
filter=none

# approximate inverse diagonal Hessian to be used at every iteration with lbfgs solver. Must have the same geometry as the model. Provide the file name. 
inverse_diagonal_hessian=none

# regularization. Possible values: -1 (no regularization), 0, 1, 2 for 0th, first and second-order Tikhonov
regularization=-1

# prior model to be used in the regularization. Assumed zero if not provided.
prior=none

# regularization total damping parameter
lambda=0

# regularization directional weights
reg_xweight=1.0
reg_zweight=1.0

# apply trace-by-trace normalization. Possible values: 0 or 1
normalize=0

# envelop inversion. Possible values: 0 or 1
envelop=0

# model parameterization. 0: lame parameters-rho-anisotropy, 1: Vp-Vs-rho-aniso, 2: IP-IS-rho-aniso, 3: log(Vs/Vs0)-log(Vp/Vs - sqrt(2))-log(rho/rho0)-aniso
model_parameterization=1

# use soft clip operator based on the model bounds defined above. 0: no soft clip (hard clipping is performed), 1: use soft clip (recommended for the inversion)
soft_clip=1

# use or not cubic B-spline model parameterization. Possible values: 0 or 1
bsplines=0

# number of B-spline nodes in every direction
bs_nx=3
bs_nz=3

# how frequently model, gradient, and residual are output
isave=10

#############################################
############ miscallenous ###################
#############################################

# verbosity logging level: 0, 1, 2, or 3
verbose=3

# index of first gpu device to use (for a multi-gpu node)
device=0

# maximum number of cpu threads to use in multi-threading. If 0, it will be calculated automatically.
nthreads=0

# All file format read/write, 0: SEPlib, 1: native binary
format=0

# datapath to store output binaries
datapath=none
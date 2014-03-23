# ica - Perform infomax independent component analysis 
#       the binary ica application. 
#
#   Master input settings file for running the ICA algorithm of Bell & Sejnowski 
#   (1996) and/or the extended-ICA algorithm of Lee, Girolami & Sejnowski (1998). 
#   via the stand-alone C++-coded binary. This is the master .sc file. 
#   Normally, do not alter this file. Matlab function binica() makes 
#   a customized copy in the pwd and runs the binary ica application from 
#   that copy. Original Matlab code (runica.m) by Scott Makeig with
#   code contributions by Tony Bell, Te-Won Lee, et al. C++ translation 
#   by Sigurd Enghoff, CNL / Salk Institute 7/98
#
#   The MATLAB binica() routine can be used to call binary ica from Matlab.
#   Usage:   >> [wts,sph] = binica(data,[runica() args]);
#
#   Contacts: {scott,jung}@sccn.ucsd.edu
#             {tony,terry}@salk.edu
#
# Required input variables:
# 
    DataFile     XXX       # Input data to decompose (native floats)
                           #  multiplexed by channel (i.e., c1, c2, ...))
    chans        31        # Number of data channels (= data columns) 
    frames       768       # Number of data points (= data rows)
#
# Output variables:
#
    WeightsOutFile  binica.wts  # Output ICA weight matrix (floats)
    SphereFile      binica.sph  # Output sphering matrix (floats)
#
# Note: input data files must be, and output files will be native floats.
#
# Processing options:
# 
    sphering     on        # Flag sphering of data (on/off)   {default: on}
    bias         on        # Perform bias adjustment (on/off) {default: on}
    extended     0         # Perform "extended-ICA" using tnah() with kurtosis
                           #  estimation every N training blocks. If N < 0,
                           #  fix number of sub-Gaussian components to -N 
                           #  {default|0: off}
    pca          0         # Decompose a principal component subspace of
#                          #  the data. Retain this many PCs. {default|0: all}
# Optional input variables:
# 
    lrate        1.0e-4    # Initial ICA learning rate (float << 1)
                           #  {default: heuristic ~5e-4}
    blocksize    0         # ICA block size (integer << datalength) 
                           #  {default: heuristic fraction of log data length}
    stop         1.0e-7    # Stop training when weight-change < this value
                           #  {default: heuristic ~0.0000001}
    maxsteps     512       # Max. number of ICA training steps {default: 128}
    posact       off       # Make each component activation net-positive 
                           # (on/off) NB: requires extra space! {default: off}
    annealstep   0.98      # Annealing factor (range (0,1]) - controls 
                           #  the speed of convergence.
    annealdeg    60        # Angledelta threshold for annealing {default: 60}
    momentum     0         # Momentum gain (range [0,1])      {default: 0}
    verbose      on        # Give ascii messages (on/off) {default: on}
#
# Optional input starting weights:
#
#   WeightsInFile  file    # Starting ICA weight matrix (nchans,ncomps)
#                          #  {default: identity or sphering matrix}
# 
# Optional output files:
# 
#  ActivationsFile data.act # Activations of each component (ncomps,points)
#  BiasFile      data.bs   # Bias weights (ncomps,1)
#  SignFile      data.sgn  # Signs designating (-1) sub- and (1) super-Gaussian 
                           #  components (ncomps,1) from 'extended' mode training
# Unused flags 
#
#   epochs       436       # Number of epochs
#   FrameWindow  20        # Number of frames per window
#   FrameStep    4         # Number of frames to step per window
#   EpochWindow  100       # Number of epochs per window
#   EpochStep    25        # Number of epochs to step per window
#   Baseline     25        # Number of data points contained in baseline

DSS MATLAB package

  Copyright (C) 2004, 2005 DSS MATLAB package team (dss@cis.hut.fi).
  Distributed by Laboratory of Computer and Information Science,
  Helsinki University of Technology. http://www.cis.hut.fi/projects/dss/.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

REQUIREMENTS

This is a MATLAB package. The command line version requires MATLAB v6.1 
and the graphical user interface v7. 

INSTALLATION

Copy archive into proper directory and unpack it

> gtar -xvzf dss_package.tar.gz

Add DSS package (directory where 'dss.m' is located) to the Matlab path.
Use absolute path definition.
   a) Use Matlab pathtool to add path.
OR
   b) Add line to the startup.m in your home directory
        addpath('/path/to/dss')
OR
   c) Add item to the p vector in ~/.matlab/pathdef.m
        '/path/to/dss:',...

If you want to use keyboard to interrupt calculation (so that it can
be continued afterwards), MEX or java function must be set to handle
keypresses.  To set up java keyboard interrupt modify 'classpath.txt'
file and add path to dss directory on a new line or set the path using
MATLABs javaclasspath command. To set up MEX keyboard interrupt a
testkeypress-function must be compiled with 'prepare.m' script. This
requires that MEX compiler environment is properly set up.


GRAPHICAL USER INTERFACE

DSS package contains a graphical user interface (GUI). See README in
gui-directory for more information. 

USAGE

DSS can be run with default parameters simply by passing the mixed
signal matrix as a single parameter. The result is the state structure
containing the details of the run. The state structure allows long
calculations to be interrupted and continued later on. The algorithm
also copies the essential results to the other output arguments.

>> [state] = denss(X)
>> [state, W, A, ...] = denss(X)

The algorithm can be interrupted by a control user interface or
pressing the SPACE-key and choosing the proper interrupt option.
After interruption, the calculation can be continued
by giving the state structure as a parameter.

>> [state, ...] = denss(state)

The parameters are given to the algorithm in a parameter structure that is
similar to the state structure. Parameters can be also given as
a cell array containing parameter and value pairs.
The initialization of the state script can be done in
a separate initialization function.
Also an existing state structure can be updated with new parameters.

% create new parameter structure, set the source dimension to 3
>> params.sdim = 3

% use symmetric approach
>> params.algorithm = 'symm'

% run algorithm with existing parameter structure
>> [state, ...] = denss(X, params)

% run algorithm by passing parameters directly
>> [state, ...] = denss(X, {'sdim', 3, 'algorithm', 'symm', ...})

% create algorithm state structure from parameter structure
>> state = dss_create_state(X, params)

% run algorithm with previously created state structure or
% continue interrupted calculation
>> [state, ...] = denss(state)

% update the algorithm state structure with new parameters
>> state = dss_create_state(state, params)

% run algorithm with existing state structure together with new parameters
>> [state, ...] = denss(state, params)

Some description of the parameters can be found from the files denss.m
and dss_init.m. DSS algorithm customization is made mostly by
providing custom functions for various algorithm operations.

See the demos in the demos directory for more information.

NOTES

MATLAB uses a clever pass-by-reference-or-value in function
calls. This means that passed variables are not copied (which would
take some CPU-time) if the function does not change them. For this
reason, one should avoid changing e.g. the state variable in the
denoising functions, even if it is not returned (no change would be
visible in the main program). Instead, make local copies of the
relevant parts in the denoising function.

If you have run the dss algorithm once and would like to run it again,
you have to clear the fields 'S' and 'W' from the state
struct. Otherwise, dss will start from the predefined S (or W).

FILES

--Main DSS functions:
denss.m                 Main script for running DSS algorithm. Creates
                        state structure and calls defl/symm implementation
dss_fastica.m           DSS with FastICA-like interface
dss_2dmask.m            2-dim mask denoising (eg. spectrogram denoising)
dss_create_state.m      Creates DSS state struct based on given parameters
dss_set_denoising.m     Initialize denoising function and parameters
estimate_mask.m         Estimate bit mask based on signal dynamics
prepare.m               Initialize environment (compile testkeypress.mex)

--Internal functions:
dss_preprocess.m        Performs data preprocessing
dss_check_adaptivity.m  Initialize alpha, beta and gamma values or functions
dss_core_defl.m         Deflation DSS core
dss_core_symm.m         Symmetric DSS core
dss_core_pca.m          PCA DSS core
message.m               Print message depending on verbosity level

--Helper functions:
randlap.m               Random laplacian noise
logplot.m               Plot data in linear & logarithmic scale
Keytest.class           Keyboard reading for interrupting DSS from keyboard

--Preprocessing functions:
pre_sphere.m            Default whitening function (basic PCA whitening)
pre_sphere_symm.m       Symmetric sphering

--Stopping functions:
default_stop.m          Default stopping criteria function

--Orthogonalization functions:
ortho_default.m         Default orthogonalization function (both defl & symm)
ortho_quasi.m           Quasiorthogonalization

--Denoising functions:

denoise_fica_gauss.m    Equals FastICA gaussian non-linearity
denoise_fica_kurtosis.m Equals FastICA kurtosis non-linearity
denoise_fica_skew.m     Equals FastICA skewness non-linearity
denoise_fica_tanh.m     Equals FastICA tanh non-linearity (default)
denoise_pow3.m          Kurtosis denoising function
denoise_tanh.m          Tanh denoising function
denoise_mask.m          Simple mask denoising function
denoise_energy.m        Energy based denoising
denoise_smooth_tanh.m   Smooth tanh denoising
denoise_filter.m        Generic signal filter (also a denoising function)
denoise_dct.m
denoise_avg.m           Ensemble average denoising using given triggers

--Alpha functions:
*none

--Beta functions:
beta_global.m           Global spectral shift for all denoising functions
beta_tanh.m             Adaptive (local) spectral shift for tanh denoising
beta_pow3.m             Adaptive (local) spectral shift for pow3
                        denoising, equals kurtosis extremisation

--Gamma functions:
gamma_179.m             179-rule adaptive gamma
gamma_predictive.m      Predictive controller adaptive gamma (defl dss)
gamma_predictive_symm.m Predictive controller adaptive gamma (symm dss)

--Reporting functions:
report_convergence.m
report_objective.m
report_print.m
report_w.m

--Directories:
gui                     DSS GUI
src                     Misc source files
test                    Unit tests
demos                   Demonstration data and scripts


CALLBACK FUNCTION SIGNATURES

Preprocessing:
  [params, Y, wM, dwM] = preprocf(params, X, dim)
    params  Function specific parameters
    X       Original data
    Y       Whitened data
    wM      Whitening matrix
    dwM     Dewhitening matrix
    dim     Optional dimension reduction

Orthogonalization:
  [params, W] = orthof(params,W)      for symmetric dss
  [params, w] = orthof(params, W, w)  for deflation dss
    params  Function specific parameters
    W       Unmixing projection matrix estimate
    w       Projection vector estimate

Stopping:
  [stop, params]  = stopf(params, state)
    state   DSS algorithm state
    params  Stopping function specific modifiable parameters
    stop    Boolean, true to stop iteration

Denoising:
  [params, s_new] = denf(params, s, state)
    params  Denoising function specific modifiable parameters
    s       Source signal estimate, matrix of row vector signals
    state   DSS algorithm state
    s_new   Denoised signal estimate

  [params, alpha] = alphaf(params, state)
  [params, beta] = betaf (params, state)
  [params, gamma] = gammaf(params, state)
    params  Coefficient function specific modifiable parameters
    state   DSS algorithm state
    alpha/beta/gamma   Coefficient variable 

Reporting:
  [report_data]   = reportf(report_data, state)
    report_data  Reporting function specific modifiable parameters
    state        DSS algorithm state


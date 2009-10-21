function [source] = sourceanalysis(cfg, data, baseline);

% SOURCEANALYSIS performs beamformer dipole analysis on EEG or MEG data
% after preprocessing and a timelocked or frequency analysis
%
% Use as either
%   [source] = sourceanalysis(cfg, freq)
%   [source] = sourceanalysis(cfg, timelock)
%
% where the data in freq or timelock should be organised in a structure
% as obtained from the FREQANALYSIS or TIMELOCKANALYSIS function. The
% configuration "cfg" is a structure containing information about
% source positions and other options.
%
% The different source reconstruction algorithms that are implemented
% are
%   cfg.method     = 'lcmv'    linear constrained minimum variance beamformer
%                    'sam'     synthetic aperture magnetometry
%                    'dics'    dynamic imaging of coherent sources
%                    'pcc'     partial cannonical correlation/coherence
%                    'mne'     minimum norm estimation
%                    'loreta'  minimum norm estimation with smoothness constraint
%                    'rv'      scan residual variance with single dipole
%                    'music'   multiple signal classification
%                    'mvl'   multivariate Laplace source localization
% The DICS and PCC methods are for frequency domain data, all other methods 
% are for time domain data.
%
% The positions of the sources can be specified as a regular 3-D
% grid that is aligned with the axes of the head coordinate system
%   cfg.grid.xgrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.ygrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.zgrid      = vector (e.g.   0:1:20) or 'auto' (default = 'auto')
%   cfg.grid.resolution = number (e.g. 1 cm) for automatic grid generation
% Alternatively the position of a few sources at locations of interest can
% be specified, for example obtained from an anatomical or functional MRI
%   cfg.grid.pos        = Nx3 matrix with position of each source
%   cfg.grid.dim        = [Nx Ny Nz] vector with dimensions in case of 3-D grid (optional)
%   cfg.grid.inside     = vector with indices of the sources inside the brain (optional)
%   cfg.grid.outside    = vector with indices of the sources outside the brain (optional)
% You can also use the PREPARE_LEADFIELD function to create a grid with
% dipole positions and with precomputed leadfields.
%
% The following strategies are supported to obtain statistics for the source parameters using
% multiple trials in the data, either directly or through a resampling-based approach
%   cfg.singletrial   = 'no' or 'yes'   construct filter from average, apply to single trials
%   cfg.rawtrial      = 'no' or 'yes'   construct filter from single trials, apply to single trials
%   cfg.jackknife     = 'no' or 'yes'   jackknife resampling of trials
%   cfg.pseudovalue   = 'no' or 'yes'   pseudovalue resampling of trials
%   cfg.bootstrap     = 'no' or 'yes'   bootstrap resampling of trials
%   cfg.numbootstrap  = number of bootstrap replications (e.g. number of original trials)
% If none of these options is specified, the average over the trials will
% be computed prior to computing the source reconstruction.
%
% To obtain statistics over the source parameters between two conditions, you
% can also use a resampling procedure that reshuffles the trials over both
% conditions. In that case, you should call the function with two datasets
% containing single trial data like
%   [source] = sourceanalysis(cfg, freqA, freqB)
%   [source] = sourceanalysis(cfg, timelockA, timelockB)
% and you should specify
%   cfg.randomization      = 'no' or 'yes'
%   cfg.permutation        = 'no' or 'yes'
%   cfg.numrandomization   = number, e.g. 500
%   cfg.numpermutation     = number, e.g. 500 or 'all'
%
% You should specify the volume conductor model with
%   cfg.hdmfile       = string, file containing the volume conduction model
% or alternatively
%   cfg.vol           = structure with volume conduction model
%
% If the sensor information is not contained in the data itself you should
% also specify the sensor information using
%   cfg.gradfile      = string, file containing the gradiometer definition
%   cfg.elecfile      = string, file containing the electrode definition
% or alternatively
%   cfg.grad          = structure with gradiometer definition
%   cfg.elec          = structure with electrode definition
%
% If you have not specified a grid with pre-computed leadfields,
% the leadfield for each grid location will be computed on the fly.
% In that case you can modify the leadfields by reducing the rank
% (i.e.  remove the weakest orientation), or by normalizing each
% column.
%   cfg.reducerank  = 'no', or number (default = 3 for EEG, 2 for MEG)
%   cfg.normalize   = 'no' or 'yes' (default = 'no')
%
% Other configuration options are
%   cfg.channel       = Nx1 cell-array with selection of channels (default = 'all'),
%                       see CHANNELSELECTION for details
%   cfg.frequency     = single number (in Hz)
%   cfg.latency       = single number in seconds, for time-frequency analysis
%   cfg.lambda        = number or empty for automatic default
%   cfg.refchan       = reference channel label (for coherence)
%   cfg.refdip        = reference dipole location (for coherence)
%   cfg.supchan       = suppressed channel label(s)
%   cfg.supdip        = suppressed dipole location(s)
%   cfg.keeptrials    = 'no' or 'yes'
%   cfg.keepleadfield = 'no' or 'yes'
%   cfg.projectnoise  = 'no' or 'yes'
%   cfg.keepfilter    = 'no' or 'yes'
%   cfg.keepcsd       = 'no' or 'yes'
%   cfg.keepmom       = 'no' or 'yes'
%   cfg.feedback      = 'no', 'text', 'textbar', 'gui' (default = 'text')
%
% See also SOURCEDESCRIPTIVES, SOURCESTATISTICS, PREPARE_LEADFIELD

% Undocumented local options:
% cfg.numcomponents
% cfg.refchannel
% cfg.trialweight        = 'equal' or 'proportional'
% cfg.powmethod          = 'lambda1' or 'trace'
%
% This function depends on PREPARE_DIPOLE_GRID which has the following options:
% cfg.grid.xgrid (default set in PREPARE_DIPOLE_GRID: cfg.grid.xgrid = 'auto'), documented
% cfg.grid.ygrid (default set in PREPARE_DIPOLE_GRID: cfg.grid.ygrid = 'auto'), documented
% cfg.grid.zgrid (default set in PREPARE_DIPOLE_GRID: cfg.grid.zgrid = 'auto'), documented
% cfg.grid.resolution, documented
% cfg.grid.pos, documented
% cfg.grid.dim, documented
% cfg.grid.inside, documented
% cfg.grid.outside, documented
% cfg.mri
% cfg.mriunits
% cfg.smooth
% cfg.sourceunits
% cfg.threshold
% cfg.symmetry
%
% This function depends on PREPARE_FREQ_MATRICES which has the following options:
% cfg.channel (default set in SOURCEANALYSIS: cfg.channel = 'all'), documented
% cfg.dicsfix
% cfg.frequency, documented
% cfg.latency, documented
% cfg.refchan, documented
%
% This function depends on PREPARE_RESAMPLED_DATA which has the following options:
% cfg.jackknife (default set in SOURCEANALYSIS: cfg.jackknife = 'no'), documented
% cfg.numbootstrap, documented
% cfg.numcondition (set in SOURCEANALYSIS: cfg.numcondition = 2)
% cfg.numpermutation (default set in SOURCEANALYSIS: cfg.numpermutation = 100), documented
% cfg.numrandomization (default set in SOURCEANALYSIS: cfg.numrandomization = 100), documented
% cfg.permutation (default set in SOURCEANALYSIS: cfg.permutation = 'no'), documented
% cfg.pseudovalue (default set in SOURCEANALYSIS: cfg.pseudovalue = 'no'), documented
% cfg.randomization (default set in SOURCEANALYSIS: cfg.randomization = 'no'), documented
%
% This function depends on PREPARE_VOL_SENS which has the following options:
% cfg.channel, (default set in SOURCEANALYSIS: cfg.channel = 'all'), documented
% cfg.elec, documented
% cfg.elecfile, documented
% cfg.grad, documented
% cfg.gradfile, documented
% cfg.hdmfile, documented
% cfg.order
% cfg.vol, documented
%
% This function depends on PREPARE_LEADFIELD which has the following options:
% cfg.feedback, (default set in SOURCEANALYSIS: cfg.feedback = 'text'), documented
% cfg.grid, documented
% cfg.lbex
% cfg.normalize (default set in SOURCEANALYSIS), documented
% cfg.previous
% cfg.reducerank (default set in SOURCEANALYSIS), documented
% cfg.sel50p
% cfg.version

% Copyright (c) 2003-2008, Robert Oostenveld, F.C. Donders Centre
%
% $Log: sourceanalysis.m,v $
% Revision 1.141  2009/10/12 14:44:07  jansch
% built in possibility (undocumented) to efficiently project single trial estimates
% through precomputed lcmv filters
%
% Revision 1.140  2009/06/04 13:37:28  marvger
% changed mvlap name to mvl to make it consistent with literature
%
% Revision 1.139  2009/05/20 16:53:30  marvger
% added support for mvlap method (tentative)
%
% Revision 1.138  2009/03/26 13:32:12  roboos
% added SAM as beamformer method
% fixed bug in the default assignment of reducerank, which for some MEG systems caused the reducerank default to be 3 instead of 2 (which is preferred)
%
% Revision 1.137  2009/03/11 11:27:34  roboos
% added a comment and removed a now obsolete channelselection call
%
% Revision 1.136  2009/03/06 13:02:59  sashae
% changed default cfg.projectnoise='yes' to 'no'
%
% Revision 1.135  2009/02/06 10:53:10  jansch
% added experimental possibility to apply prewhitening. fixed incorrect
% conjugate transposition and replaced by explicit transpose()
%
% Revision 1.134  2009/02/05 10:22:07  roboos
% better support for single-trial data, thanks to Vladimir
%
% Revision 1.133  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.132  2009/01/19 12:29:51  roboos
% preallocate the fake covariance matrix in case it is absent in the data
% ensure that length(size(cov))==2 in case of Ntrials==1 (i.e. when keeptrials=yes and only one trial in the data)
%
% Revision 1.131  2008/11/21 13:21:35  sashae
% added call to checkconfig at start and end of fucntion
%
% Revision 1.130  2008/10/02 15:32:21  sashae
% replaced call to createsubcfg with checkconfig
%
% Revision 1.129  2008/09/26 12:42:18  sashae
% checkconfig: checks if the input cfg is valid for this function
%
% Revision 1.128  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.127  2008/07/15 19:56:44  roboos
% moved cfg details for dipole grid to subcfg (cfg.grid)subcfg (cfg.grid.xxx)
%
% Revision 1.126  2008/07/07 12:37:42  roboos
% fixed bug in channelordering in case sens.label and data.label were inconsistent
%
% Revision 1.125  2008/04/10 08:03:11  roboos
% renamed the fieldtrip/private/prepare_vol_sens function into prepare_headmodel
%
% Revision 1.124  2007/05/15 07:01:29  roboos
% updated help
%
% Revision 1.123  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.122  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.121  2007/03/27 11:05:19  ingnie
% changed call to fixdimord in call to checkinput
%
% Revision 1.120  2007/03/05 15:31:37  roboos
% added empty defaults for refchan/supchan in case of method=pcc
%
% Revision 1.119  2006/10/19 15:19:14  roboos
% switched to using beamformer_lcmv and beamformer_dics for the specific methods
%
% Revision 1.118  2006/10/12 13:04:34  roboos
% updated documentation
%
% Revision 1.117  2006/10/04 08:10:51  roboos
% changed default for cfg.reducerank, new default is 2 for MEG and 3 for EEG
%
% Revision 1.116  2006/10/04 07:10:07  roboos
% updated documentation
%
% Revision 1.115  2006/10/02 13:01:21  roboos
% also recognize rpttap as part of the dimord
%
% Revision 1.114  2006/08/16 08:23:36  roboos
% updated documentation
%
% Revision 1.113  2006/07/05 10:34:52  roboos
% minor change in documentation
%
% Revision 1.112  2006/07/05 10:32:14  roboos
% minor change in documentation
%
% Revision 1.111  2006/07/05 10:26:11  roboos
% updated documentation, the cfg.xgrid/ygrid/zgrid options have moved to cfg.grid substructure
% removed default 'auto' values for xgrid/ygrid/zgrid (were similar to the default values in prepare_dipole_grid)
%
% Revision 1.110  2006/07/04 17:06:23  ingnie
% updated documentation
%
% Revision 1.109  2006/07/04 16:04:50  roboos
% renamed option 'jacknife' into 'jackknife' for consistency, maintain backward compatibility with cfgs and old data
%
% Revision 1.108  2006/06/26 09:26:03  ingnie
% fixed bug (forgotten "end")
%
% Revision 1.107  2006/06/22 12:27:41  roboos
% optionally pass the covariance as input argument to MUSIC
%
% Revision 1.106  2006/06/20 16:25:56  ingnie
% updated documentation
%
% Revision 1.105  2006/06/14 11:55:07  roboos
% switched from construct_optarg() subfunction to the new createsubcfg function
%
% Revision 1.104  2006/06/13 14:48:09  ingnie
% updated documentation
%
% Revision 1.103  2006/05/23 10:15:31  roboos
% Changed the initial construction of the dipole grid, either with or
% without precomputed leadfields (there was a case in which the automatic
% grid would not be tight). Improved documentation around the disabled
% section for singletrial=yes. Some other small changes.
%
% Revision 1.102  2006/05/10 08:14:09  roboos
% removed support for parallelization on cluster, only precompute leadfields
% if not already present, implemented subfunction construct_optarg() for
% clean and consistent construction of key-value arguments for the low-level
% inverse methods, allow additional arguments to be passed transparently
% to the low-level inverse methods by means of cfg.xxx.key=value, where
% xxx=pcc/lcmv/dics/rv/music/loreta/mne
%
% Revision 1.101  2006/05/03 15:10:11  roboos
% only define submethod in the case of method=dics
%
% Revision 1.100  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.99  2006/04/12 08:38:01  ingnie
% updated documentation
%
% Revision 1.98  2006/04/10 16:33:46  ingnie
% updated documentation
%
% Revision 1.97  2006/04/06 16:17:10  ingnie
% updated documentation
%
% Revision 1.96  2006/03/22 07:53:06  roboos
% small change in documentation
%
% Revision 1.95  2006/03/09 08:19:06  roboos
% apply fixdimord to the baseline condition data if present
%
% Revision 1.94  2006/02/24 16:41:39  roboos
% changed foi and toi into freq and time for frequency data
%
% Revision 1.93  2006/02/23 10:28:16  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.92  2006/02/07 22:21:01  roboos
% changed the xgrid/ygrid/zgrid and dim in the output source structure (from cfg to grid), they can be different from the cfg since prepare_dipole_grid will make the box tight
%
% Revision 1.91  2006/02/07 20:08:23  roboos
% changed all occurences of a dimord with chancmb (was previous sgncmb) into chan
%
% Revision 1.90  2006/02/01 12:26:00  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.89  2006/01/24 20:51:57  roboos
% fixed obvious typo
%
% Revision 1.88  2006/01/24 20:02:50  roboos
% renamed dics_cohrefdip and cohrefchan into cfg.method=dics (backwards compatible)
% added extra local variable "submethod" that indicates the details for dics
%
% Revision 1.87  2006/01/24 14:28:33  roboos
% small change in documentation
%
% Revision 1.86  2005/11/24 16:25:16  roboos
% added average or single trial ERF to the input of PCC beamformer in timedomain
%
% Revision 1.85  2005/11/16 09:07:03  roboos
% added option cfg.powmethod to determine whether lambda1 or trace is used in beamformer
%
% Revision 1.84  2005/11/08 11:07:13  roboos
% added support for normalize and reducerank, by passing these options on to the beamformer subfunction
%
% Revision 1.83  2005/09/29 00:45:00  roboos
% added multiple signal classification as reconstruction technique for timelocked data
%
% Revision 1.82  2005/10/25 08:42:56  roboos
% changed the detection of either frequency or timelock data using flags (isfreq, iscomp, istimelock) instead of by renaming the input data object to freq/timelock
% renamed the freq/timelock variable throughout the code into "data"
% added support for freq/comp data to the non-beamforming timelocked source reconstruction methods
%
% Revision 1.81  2005/10/14 15:47:03  roboos
% add an identity covariance matrix to timelocked data if not present
% some cosmetic changes
%
% Revision 1.80  2005/10/05 11:14:11  roboos
% updated documentation
% added support for LORETA (still experimental)
% removed all code that referred to unimplemented parallelization for some of the reconstruction methods
%
% Revision 1.79  2005/09/05 06:36:56  jansch
% keep cumtapcnt of freq-input in output
%
% Revision 1.78  2005/08/16 13:15:55  jansch
% *** empty log message ***
%
% Revision 1.77  2005/08/16 12:43:03  jansch
% included possibility to have 'rpttap_sgncmb_frq' as input
%
% Revision 1.76  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.75  2005/06/17 09:01:13  roboos
% moved implementation of keepmom and keepcsd over to beamformer subfunction
%
% Revision 1.74  2005/05/17 17:56:38  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
% implemented support for projecting the single trial fft matrix through the filters of beamformer_pcc
%
% Revision 1.73  2005/04/25 16:05:11  roboos
% includedupdated options for cortex segmentation based grid in the documentation (not in main documentation)
%
% Revision 1.72  2005/03/03 13:57:29  roboos
% implemented pseudovalue resampling, i.e. a combination of jackknife and the complete average
%
% Revision 1.71  2005/03/03 10:57:43  roboos
% Renamed the source analysis methods, the submethods for dics are now
% recognized by their additional required fields (refchan/refdip). Added
% support for experimental pcc method (not yet implemented in a clean way),
% pcc is now a re3cognised submethod of the standard beamformer function.
% Included mne method in the documentation.
%
% Revision 1.70  2005/02/14 21:53:14  roboos
% changed assignment of empty struct to dip into 'clear dip' (for mentat-style parallelization)
%
% Revision 1.69  2005/01/26 08:38:50  roboos
% extended and fixed the experimental partcanoncorr code
%
% Revision 1.68  2005/01/11 12:43:03  roboos
% added temporary hack to support partcanoncorr as a beamformer method
%
% Revision 1.67  2005/01/10 17:26:46  roboos
% added support for 2nd method of parallelizing the code, separated
% two methods by the name of the backend toolbox (beowulf or mentat)
%
% Revision 1.66  2004/12/08 18:00:13  roboos
% implemented consistent method of selecting a subset of channels for
% forward and inverse computations using cfg.channel and updated the
% ducumentation
%
% Revision 1.65  2004/10/21 09:01:34  roboos
% added numcondition=2 for averaging two conditions while resampling
%
% Revision 1.64  2004/09/28 08:01:03  roboos
% added residual variance scanning with single dipole (cfg.method=rv)
% added 'yes' as allowed value for cfg.parallel
%
% Revision 1.63  2004/09/21 14:30:07  roboos
% fixed bug in DICS/power for two conditions (Nbaseline~=Ntrials)
%
% Revision 1.62  2004/09/21 13:22:15  roboos
% in the previous revision (1.61) initial support for MNE was added but accidentally not logged in CVS
% this revision only adds a single comment line with the undocumented method=mne option
%
% Revision 1.61  2004/09/21 13:15:32  roboos
% *** empty log message ***
%
% Revision 1.60  2004/09/21 12:49:00  roboos
% fixed a bug in averaging of 2 input conditions using prepare_resampled_data (using cfg2)
%
% Revision 1.59  2004/09/14 16:28:49  roboos
% incorporated the changes on the "exp-resampling-parallel" branch
% new implementation of resampling/randomization
% new implementation of parallellization
%
% Revision 1.58.4.5  2004/09/08 16:36:00  roboos
% fixed many bugs that stopped execution through rigorous testing
% correctness of the output still has to be evaluated
% non-proportional trial weighing is not implemented yet
%
% Revision 1.58.4.4  2004/09/08 12:25:13  roboos
% auto-indented code using matlab editor, fixed one incorrect "end"
%
% Revision 1.58.4.3  2004/09/08 09:35:49  roboos
% implemented prepare_resampled_data and new evalwulf for DICS, some related changes to LCMV
%
% Revision 1.58.4.2  2004/09/06 16:32:28  roboos
% implemented new prepare_resampled_data for a few resampling strategies and lcmv only, untested
%
% Revision 1.58.4.1  2004/09/06 16:10:26  roboos
% disabled parallelization over grid (still in code)
% moved documentation for paralellization to hidden location
% fixed bug in jackknife+lcmv
%
% Revision 1.58  2004/08/16 12:14:33  roboos
% fixed small bug that was introduced in the previous update
%
% Revision 1.57  2004/08/16 10:53:09  roboos
% fixed bug for lcmv together with randomization (trial assignment was oubviously wrong)
%
% Revision 1.56  2004/08/16 10:01:23  roboos
% updated help, fixed bug in output assignment with permutation, changed handling of keepleadfield and precomputed leadfields
%
% Revision 1.55  2004/08/06 12:23:04  roboos
% fixed small bug in default setting for cfg.permutation
%
% Revision 1.54  2004/08/05 07:14:25  roboos
% implemented permutation analogous to randomization, updated and restructured help
%
% Revision 1.53  2004/08/05 06:37:29  roboos
% extra removal of leadfield in output source and cfg, removed old revision comments
%
% Revision 1.52  2004/08/04 11:53:03  roboos
% updated the help
%
% Revision 1.51  2004/06/09 10:07:19  roberto
% added version information and input data configuration to output cfg
%
% Revision 1.50  2004/03/01 11:01:55  roberto
% added precomputed filters to grid and using these, implemented singletrial for dics
%
% Revision 1.49  2004/02/05 12:19:36  roberto
% moved the code that was common with dipolefitting into separate prepare_grid and prepare_sens_vol functions
% also renamed internal subfunction freq_to_matrices to prepare_freq_matrices and moved into separate function
%
% Revision 1.48  2004/02/02 10:48:11  roberto
% fixed two bugs in randomization for lcmv & timelocked data
%
% Revision 1.47  2004/01/26 09:04:54  roberto
% fixed serious bug in selection of spheres for multi-sphere volume model
% selection was based on wrong indices, causing reference magnetometer spheres
% and gradiometer spheres to be mixed up
%
% Revision 1.46  2004/01/22 10:50:55  roberto
% fixed bug in cfg.keepmom
% improved handling of pre-computation of leadfield without scanning
%
% Revision 1.45  2004/01/21 13:07:52  roberto
% added a fprintf-line to indicate the progress in randomization (which takes quite some time)
%
% Revision 1.44  2004/01/20 09:09:45  roberto
% added computation of leadfield in case not provided and multiple trials have to
% be scanned
%
% Revision 1.43  2004/01/14 14:29:43  roberto
% changed handling of sgncmbindx etc inside freq_to_matrices subfunction
% mainly cosmetical changes, but also fixed bug in refindx for pow_refchan
%
% Revision 1.42  2004/01/14 08:50:05  roberto
% changed the name of output variables for randomization
%
% Revision 1.41  2004/01/08 09:25:35  roberto
% canged cfg.label into freq.sgn inside freq_to_matrices subfunction to fix bug
%
% Revision 1.40  2004/01/08 09:14:49  roberto
% added option keepmom to save memory, default is 'yes'

fieldtripdefs

% set a timer to determine how long the sourceanalysis takes in total
tic;

cfg = checkconfig(cfg, 'trackconfig', 'on');

% check if the input data is valid for this function
data = checkdata(data, 'datatype', {'timelock', 'freq', 'comp'}, 'feedback', 'yes');
if nargin>2
  baseline = checkdata(baseline, 'datatype', {'timelock', 'freq', 'comp'}, 'feedback', 'yes');
end

% check if the input cfg is valid for this function
cfg = checkconfig(cfg, 'renamed',     {'jacknife',   'jackknife'});
cfg = checkconfig(cfg, 'renamed',     {'refchannel', 'refchan'});
cfg = checkconfig(cfg, 'renamedval',  {'method', 'power',           'dics'});
cfg = checkconfig(cfg, 'renamedval',  {'method', 'coh_refchan',     'dics'});
cfg = checkconfig(cfg, 'renamedval',  {'method', 'coh_refdip',      'dics'});
cfg = checkconfig(cfg, 'renamedval',  {'method', 'dics_cohrefchan', 'dics'});
cfg = checkconfig(cfg, 'renamedval',  {'method', 'dics_cohrefdip',  'dics'});
cfg = checkconfig(cfg, 'forbidden',   {'parallel'});

% determine the type of input data
if isfield(data, 'freq')
  isfreq     = 1;
  iscomp     = 0;
  istimelock = 0;
elseif isfield(data, 'topo')
  isfreq     = 0;
  iscomp     = 1;
  istimelock = 0;
elseif isfield(data, 'time')
  iscomp     = 0;
  isfreq     = 0;
  istimelock = 1;
else
  error('input data is not recognized');
end

% set the defaults
if ~isfield(cfg, 'method') && istimelock, cfg.method = 'lcmv';      end
if ~isfield(cfg, 'method') && isfreq,     cfg.method = 'dics';      end
if ~isfield(cfg, 'keeptrials')        cfg.keeptrials  = 'no';     end
if ~isfield(cfg, 'keepfilter')        cfg.keepfilter = 'no';      end
if ~isfield(cfg, 'keepleadfield')     cfg.keepleadfield = 'no';   end
if ~isfield(cfg, 'keepcsd')           cfg.keepcsd     = 'no';     end
if ~isfield(cfg, 'keepmom')           cfg.keepmom     = 'yes';    end
if ~isfield(cfg, 'projectnoise')      cfg.projectnoise = 'no';    end
if ~isfield(cfg, 'trialweight')       cfg.trialweight = 'equal';  end
if ~isfield(cfg, 'jackknife'),        cfg.jackknife    = 'no';    end
if ~isfield(cfg, 'pseudovalue'),      cfg.pseudovalue = 'no';     end
if ~isfield(cfg, 'bootstrap'),        cfg.bootstrap   = 'no';     end
if ~isfield(cfg, 'singletrial'),      cfg.singletrial = 'no';     end
if ~isfield(cfg, 'rawtrial'),         cfg.rawtrial    = 'no';     end
if ~isfield(cfg, 'randomization'),    cfg.randomization = 'no';   end
if ~isfield(cfg, 'numrandomization'), cfg.numrandomization = 100; end
if ~isfield(cfg, 'permutation'),      cfg.permutation = 'no';     end
if ~isfield(cfg, 'numpermutation'),   cfg.numpermutation = 100;   end
if ~isfield(cfg, 'wakewulf'),         cfg.wakewulf    = 'yes';    end
if ~isfield(cfg, 'killwulf'),         cfg.killwulf    = 'yes';    end
if ~isfield(cfg, 'feedback'),         cfg.feedback    = 'text';   end
if ~isfield(cfg, 'supdip'),           cfg.supdip = [];            end
if ~isfield(cfg, 'lambda'),           cfg.lambda = [];            end
if ~isfield(cfg, 'powmethod'),        cfg.powmethod = [];         end
if ~isfield(cfg, 'channel'),          cfg.channel = 'all';        end
if ~isfield(cfg, 'normalize'),        cfg.normalize = 'no';       end
if ~isfield(cfg, 'prewhiten'),        cfg.prewhiten = 'no';       end
% if ~isfield(cfg, 'reducerank'),     cfg.reducerank = 'no';      end  % the default for this depends on EEG/MEG and is set below

% put the low-level options pertaining to the source reconstruction method in their own field
% put the low-level options pertaining to the dipole grid in their own field
cfg = checkconfig(cfg, 'createsubcfg',  {cfg.method, 'grid'});

convertfreq = 0;
convertcomp = 0;
if ~istimelock && (strcmp(cfg.method, 'mne') || strcmp(cfg.method, 'loreta') || strcmp(cfg.method, 'rv') || strcmp(cfg.method, 'music'))
  % these timelock methods are also supported for frequency or component data
  if isfreq
    [data, cfg] = freq2timelock(cfg, data);
    convertfreq = 1;  % flag indicating that the data was converted
  elseif iscomp
    [data, cfg] = comp2timelock(cfg, data);
    convertcomp = 1;  % flag indicating that the data was converted
  end
  istimelock = 1;     % from now on the data can be treated as timelocked
  isfreq     = 0;
  iscomp     = 0;
end

% select only those channels that are present in the data
cfg.channel = channelselection(cfg.channel, data.label); 

if nargin>2 && (strcmp(cfg.randomization, 'no') && strcmp(cfg.permutation, 'no') && strcmp(cfg.prewhiten, 'no'))
  error('input of two conditions only makes sense if you want to randomize or permute, or if you want to prewhiten');
elseif nargin<3 && (strcmp(cfg.randomization, 'yes') || strcmp(cfg.permutation, 'yes'))
  error('randomization or permutation requires that you give two conditions as input');
end

if isfield(cfg, 'latency') && istimelock
  error('specification of cfg.latency is only required for time-frequency data');
end

if sum([strcmp(cfg.jackknife, 'yes'), strcmp(cfg.bootstrap, 'yes'), strcmp(cfg.pseudovalue, 'yes'), strcmp(cfg.singletrial, 'yes'), strcmp(cfg.rawtrial, 'yes'), strcmp(cfg.randomization, 'yes'), strcmp(cfg.permutation, 'yes')])>1
  error('jackknife, bootstrap, pseudovalue, singletrial, rawtrial, randomization and permutation are mutually exclusive');
end

if isfreq
  if ~strcmp(data.dimord, 'chan_freq')          && ...
      ~strcmp(data.dimord, 'chan_freq_time')     && ...
      ~strcmp(data.dimord, 'rpt_chan_freq')      && ...
      ~strcmp(data.dimord, 'rpt_chan_freq_time') && ...
      ~strcmp(data.dimord, 'rpttap_chan_freq')   && ...
      ~strcmp(data.dimord, 'rpttap_chan_freq_time')
    error('dimord of input frequency data is not recognized');
  end
end

% collect and preprocess the electrodes/gradiometer and head model
[vol, sens, cfg] = prepare_headmodel(cfg, data);

% It might be that the number of channels in the data, the number of
% channels in the electrode/gradiometer definition and the number of
% channels in the multisphere volume conduction model are different. 
% Hence a subset of the data channels will be used.
Nchans = length(cfg.channel);

% set the default for reducing the rank of the leadfields, this is an
% option to the specific method and will be passed on to the low-level
% function
if ~isfield(cfg.(cfg.method), 'reducerank')
  if senstype(sens, 'meg')
    cfg.(cfg.method).reducerank = 2;
  else
    cfg.(cfg.method).reducerank = 3;
  end
end

if strcmp(cfg.keepleadfield, 'yes') && (~isfield(cfg, 'grid') || ~isfield(cfg.grid, 'leadfield'))
  % precompute the leadfields upon the users request
  fprintf('precomputing leadfields\n');
  [grid, cfg] = prepare_leadfield(cfg, data);
elseif (strcmp(cfg.permutation,   'yes') || ...
    strcmp(cfg.randomization, 'yes') || ...
    strcmp(cfg.bootstrap,     'yes') || ...
    strcmp(cfg.jackknife,      'yes') || ...
    strcmp(cfg.pseudovalue,   'yes') || ...
    strcmp(cfg.singletrial,   'yes') || ...
    strcmp(cfg.rawtrial,      'yes')) && (~isfield(cfg, 'grid') || ~isfield(cfg.grid, 'leadfield'))
  % also precompute the leadfields if multiple trials have to be processed
  fprintf('precomputing leadfields for efficient handling of multiple trials\n');
  [grid, cfg] = prepare_leadfield(cfg, data);
else
  % only prepare the grid positions, the leadfield will be computed on the fly if not present
  [grid, cfg] = prepare_dipole_grid(cfg, vol, sens);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do frequency domain source reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfreq && any(strcmp(cfg.method, {'dics', 'pcc'}))

  if strcmp(cfg.method, 'pcc')
    % HACK: requires some extra defaults
    if ~isfield(cfg, 'refdip'), cfg.refdip = []; end
    if ~isfield(cfg, 'supdip'), cfg.supdip = []; end
    if ~isfield(cfg, 'refchan'), cfg.refchan = []; end
    if ~isfield(cfg, 'supchan'), cfg.supchan = []; end
    cfg.refchan = channelselection(cfg.refchan, data.label);
    cfg.supchan = channelselection(cfg.supchan, data.label);

    % HACK: use some experimental code
    if nargin>2 && strcmp(cfg.prewhiten, 'no'),
      error('not supported')
    end
    tmpcfg         = cfg;
    tmpcfg.refchan = ''; % prepare_freq_matrices should not know explicitly about the refchan
    tmpcfg.channel = cfg.channel(:)';
    if isfield(cfg, 'refchan')
      % add the refchan implicitely
      tmpcfg.channel = [tmpcfg.channel cfg.refchan(:)'];
    end
    if isfield(cfg, 'supchan')
      % add the supchan implicitely
      tmpcfg.channel = [tmpcfg.channel cfg.supchan(:)'];
    end

    % select the data in the channels and the frequency of interest
    [Cf, Cr, Pr, Ntrials, tmpcfg] = prepare_freq_matrices(tmpcfg, data);
    if strcmp(cfg.prewhiten, 'yes'), 
      [Cfb, Crb, Prb, Ntrialsb, tmpcfgb] = prepare_freq_matrices(tmpcfg, baseline);
      Cf      = prewhitening_filter2(squeeze(mean(Cf,1)), squeeze(mean(Cfb,1)));
      Ntrials = 1;
    end

    if isfield(cfg, 'refchan') && ~isempty(cfg.refchan)
      [dum, refchanindx] = match_str(cfg.refchan, tmpcfg.channel);
    else
      refchanindx = [];
    end
    if isfield(cfg, 'supchan') && ~isempty(cfg.supchan)
      [dum, supchanindx] = match_str(cfg.supchan,tmpcfg.channel);
    else
      supchanindx = [];
    end
    Nchans = length(tmpcfg.channel); % update the number of channels

    % if the input data has a complete fourier spectrum, project it through the filters
    % FIXME it was incorrect , since the 
    % ' leads to a conjugate transposition check this in beamformer_pcc
    if isfield(data, 'fourierspctrm')
      [dum, datchanindx] = match_str(tmpcfg.channel, data.label);
      fbin = nearest(data.freq, cfg.frequency);
      if strcmp(data.dimord, 'chan_freq')
        avg = data.fourierspctrm(datchanindx, fbin);
      elseif strcmp(data.dimord, 'rpt_chan_freq') || strcmp(data.dimord, 'rpttap_chan_freq'),
        %avg = data.fourierspctrm(:, datchanindx, fbin)';
        avg = transpose(data.fourierspctrm(:, datchanindx, fbin));
      elseif strcmp(data.dimord, 'chan_freq_time')
        tbin = nearest(data.time, cfg.latency);
        avg = data.fourierspctrm(datchanindx, fbin, tbin);
      elseif strcmp(data.dimord, 'rpt_chan_freq_time') || strcmp(data.dimord, 'rpttap_chan_freq_time'),
        tbin = nearest(data.time, cfg.latency);
        %avg = data.fourierspctrm(:, datchanindx, fbin, tbin)';
        avg = transpose(data.fourierspctrm(:, datchanindx, fbin, tbin));
      end
    else
      avg = [];
    end

  else
    % HACK: use the default code
    % convert the input data, so that Cf, Cr and Pr contain either the average over all trials (if Ntrials==1)
    % or the individual cross-spectral-densities/powers of each individual trial (if Ntrials>1)
    [Cf, Cr, Pr, Ntrials, cfg] = prepare_freq_matrices(cfg, data);
  end

  if strcmp(cfg.method, 'dics')
    % assign a descriptive name to each of the dics sub-methods, the default is power only
    if strcmp(cfg.method, 'dics') && isfield(cfg, 'refdip') && ~isempty(cfg.refdip);
      submethod = 'dics_refdip';
    elseif strcmp(cfg.method, 'dics') && isfield(cfg, 'refchan') && ~isempty(cfg.refchan);
      submethod = 'dics_refchan';
    else
      submethod = 'dics_power';
    end
  end

  % fill these with NaNs, so that I dont have to treat them separately
  if isempty(Cr), Cr = nan*zeros(Ntrials, Nchans, 1); end
  if isempty(Pr), Pr = nan*zeros(Ntrials, 1, 1); end

  if nargin>2
    % repeat the conversion for the baseline condition
    [bCf, bCr, bPr, Nbaseline, cfg] = prepare_freq_matrices(cfg, baseline);
    % fill these with NaNs, so that I dont have to treat them separately
    if isempty(bCr), bCr = nan*zeros(Nbaseline, Nchans, 1); end
    if isempty(bPr), bPr = nan*zeros(Nbaseline, 1, 1); end
    % rename the active condition for convenience
    aCf = Cf;
    aCr = Cr;
    aPr = Pr;
    % this is required for averaging 2 conditions using prepare_resampled_data
    cfg2 = [];
    cfg2.numcondition = 2;
    % this is required for randomizing/permuting 2 conditions using prepare_resampled_data
    cfg.numcondition = 2;
  end

  % prepare the resampling of the trials, or average the data if multiple trials are present and no resampling is neccessary
  if (Ntrials<=1) && (strcmp(cfg.jackknife, 'yes') || strcmp(cfg.bootstrap, 'yes') || strcmp(cfg.pseudovalue, 'yes') || strcmp(cfg.singletrial, 'yes') || strcmp(cfg.rawtrial, 'yes') || strcmp(cfg.randomization, 'yes') || strcmp(cfg.permutation, 'yes'))
    error('multiple trials required in the data\n');

  elseif strcmp(cfg.permutation, 'yes')
    % compute the cross-spectral density matrix without resampling
    [dum, avg_aCf, avg_aCr, avg_aPr, avg_bCf, avg_bCr, avg_bPr] = prepare_resampled_data(cfg2 , aCf, aCr, aPr, bCf, bCr, bPr);
    % compute the cross-spectral density matrix with random permutation
    [dum, rnd_aCf, rnd_aCr, rnd_aPr, rnd_bCf, rnd_bCr, rnd_bPr] = prepare_resampled_data(cfg, aCf, aCr, aPr, bCf, bCr, bPr);
    % concatenate the different resamplings
    Cf = cat(1, reshape(avg_aCf, [1 Nchans Nchans]), reshape(avg_bCf, [1 Nchans Nchans]), rnd_aCf, rnd_bCf);
    Cr = cat(1, reshape(avg_aCr, [1 Nchans 1     ]), reshape(avg_bCr, [1 Nchans 1     ]), rnd_aCr, rnd_bCr);
    Pr = cat(1, reshape(avg_aPr, [1 1      1     ]), reshape(avg_bPr, [1 1      1     ]), rnd_aPr, rnd_bPr);
    % clear temporary working copies
    clear avg_aCf avg_aCr avg_aPr avg_bCf avg_bCr avg_bPr
    clear rnd_aCf rnd_aCr rnd_aPr rnd_bCf rnd_bCr rnd_bPr
    % the order of the resamplings should be [avgA avgB rndA rndB rndA rndB rndA rndB ....]
    Nrepetitions = 2*cfg.numpermutation + 2;
    order = [1 2 3:2:Nrepetitions 4:2:Nrepetitions];
    Cf = Cf(order,:,:);
    Cr = Cr(order,:,:);
    Pr = Pr(order,:,:);

  elseif strcmp(cfg.randomization, 'yes')
    % compute the cross-spectral density matrix without resampling
    [dum, avg_aCf, avg_aCr, avg_aPr, avg_bCf, avg_bCr, avg_bPr] = prepare_resampled_data(cfg2 , aCf, aCr, aPr, bCf, bCr, bPr);
    % compute the cross-spectral density matrix with random resampling
    [dum, rnd_aCf, rnd_aCr, rnd_aPr, rnd_bCf, rnd_bCr, rnd_bPr] = prepare_resampled_data(cfg, aCf, aCr, aPr, bCf, bCr, bPr);
    % concatenate the different resamplings
    Cf = cat(1, reshape(avg_aCf, [1 Nchans Nchans]), reshape(avg_bCf, [1 Nchans Nchans]), rnd_aCf, rnd_bCf);
    Cr = cat(1, reshape(avg_aCr, [1 Nchans 1     ]), reshape(avg_bCr, [1 Nchans 1     ]), rnd_aCr, rnd_bCr);
    Pr = cat(1, reshape(avg_aPr, [1 1      1     ]), reshape(avg_bPr, [1 1      1     ]), rnd_aPr, rnd_bPr);
    % clear temporary working copies
    clear avg_aCf avg_aCr avg_aPr avg_bCf avg_bCr avg_bPr
    clear rnd_aCf rnd_aCr rnd_aPr rnd_bCf rnd_bCr rnd_bPr
    % the order of the resamplings should be [avgA avgB rndA rndB rndA rndB rndA rndB ....]
    Nrepetitions = 2*cfg.numrandomization + 2;
    order = [1 2 3:2:Nrepetitions 4:2:Nrepetitions];
    Cf = Cf(order,:,:);
    Cr = Cr(order,:,:);
    Pr = Pr(order,:,:);

  elseif strcmp(cfg.jackknife, 'yes')
    % compute the cross-spectral density matrix with jackknife resampling
    [cfg, Cf, Cr, Pr] = prepare_resampled_data(cfg, Cf, Cr, Pr);
    Nrepetitions = Ntrials;

  elseif strcmp(cfg.bootstrap, 'yes')
    % compute the cross-spectral density matrix with bootstrap resampling
    [cfg, Cf, Cr, Pr] = prepare_resampled_data(cfg, Cf, Cr, Pr);
    Nrepetitions = cfg.numbootstrap;

  elseif strcmp(cfg.pseudovalue, 'yes')
    % compute the cross-spectral density matrix with pseudovalue resampling
    [cfg, Cf, Cr, Pr] = prepare_resampled_data(cfg, Cf, Cr, Pr);
    Nrepetitions = Ntrials+1;

  elseif strcmp(cfg.singletrial, 'yes')
    % The idea is that beamformer uses the average covariance to construct the
    % filter and applies it to the single trial covariance/csd. The problem
    % is that beamformer will use the averaged covariance/csd to estimate the
    % power and not the single trial covariance/csd
    error('this option contains a bug, and is therefore not supported at the moment');
    Cf = Cf; % FIXME, should be averaged and repeated for each trial
    Cr = Cr; % FIXME, should be averaged and repeated for each trial
    Pr = Pr; % FIXME, should be averaged and repeated for each trial
    Nrepetitions = Ntrials;

  elseif strcmp(cfg.rawtrial, 'yes')
    % keep all the individual trials, do not average them
    Cf = Cf;
    Cr = Cr;
    Pr = Pr;
    Nrepetitions = Ntrials;

  elseif Ntrials>1
    % compute the average from the individual trials
    Cf = reshape(sum(Cf, 1) / Ntrials, [Nchans Nchans]);
    Cr = reshape(sum(Cr, 1) / Ntrials, [Nchans 1]);
    Pr = reshape(sum(Pr, 1) / Ntrials, [1 1]);
    Nrepetitions = 1;

  elseif Ntrials==1
    % no rearrangement of trials is neccesary, the data already represents the average
    Cf = Cf;
    Cr = Cr;
    Pr = Pr;
    Nrepetitions = 1;
  end

  % reshape so that it also looks like one trial (out of many)
  if Nrepetitions==1
    Cf  = reshape(Cf , [1 Nchans Nchans]);
    Cr  = reshape(Cr , [1 Nchans 1]);
    Pr  = reshape(Pr , [1 1 1]);
  end
  
  % get the relevant low level options from the cfg and convert into key-value pairs
  optarg = cfg2keyval(getfield(cfg, cfg.method));
  
  for i=1:Nrepetitions
    fprintf('scanning repetition %d\n', i);
    if     strcmp(cfg.method, 'dics') && strcmp(submethod, 'dics_power')
      dip(i) = beamformer_dics(grid, sens, vol, [],  squeeze(Cf(i,:,:)), optarg{:});
    elseif strcmp(cfg.method, 'dics') && strcmp(submethod, 'dics_refchan')
      dip(i) = beamformer_dics(grid, sens, vol, [],  squeeze(Cf(i,:,:)), optarg{:}, 'Cr', Cr(i,:), 'Pr', Pr(i));
    elseif strcmp(cfg.method, 'dics') && strcmp(submethod, 'dics_refdip')
      dip(i) = beamformer_dics(grid, sens, vol, [],  squeeze(Cf(i,:,:)), optarg{:}, 'refdip', cfg.refdip);
    elseif strcmp(cfg.method, 'pcc')
      dip(i) = beamformer_pcc(grid, sens, vol, avg, squeeze(Cf(i,:,:)), optarg{:}, 'refdip', cfg.refdip, 'refchan', refchanindx, 'supdip', cfg.supdip, 'supchan', supchanindx);
    else
      error(sprintf('method ''%s'' is unsupported for source reconstruction in the frequency domain', cfg.method));
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % do time domain source reconstruction
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif istimelock && any(strcmp(cfg.method, {'lcmv', 'sam', 'mne', 'loreta', 'rv', 'music', 'pcc', 'mvl'}))

  % determine the size of the data
  Nsamples = size(data.avg,2);
  Nchans   = length(data.label);
  if isfield(data, 'cov') && length(size(data.cov))==3
    Ntrials = size(data.cov,1);
  elseif isfield(data, 'trial') && length(size(data.trial))==3
    Ntrials = size(data.trial,1);
  else
    Ntrials = 1;
  end

  if isfield(data, 'cov')
    % use the estimated data covariance matrix
    hascovariance = 1;
  else
    % add a identity covariance matrix, this simplifies the handling of the different source reconstruction methods
    % since the covariance is only used by some reconstruction methods and might not allways be present in the data
    if Ntrials==1
      data.cov = eye(Nchans);
    else
      data.cov = zeros(Ntrials,Nchans,Nchans);
      for i=1:Ntrials
        data.cov(i,:,:) = eye(Nchans);
      end
    end
    hascovariance = 0;
  end

  if strcmp(cfg.method, 'pcc')
    % HACK: requires some extra defaults
    if ~isfield(cfg, 'refdip'), cfg.refdip = []; end
    if ~isfield(cfg, 'supdip'), cfg.supdip = []; end

    % HACK: experimental code
    if nargin>2
      error('not supported')
    end
    tmpcfg = [];
    tmpcfg.channel = cfg.channel(:)';
    if isfield(cfg, 'refchan')
      tmpcfg.channel = [tmpcfg.channel cfg.refchan(:)'];
    end
    if isfield(cfg, 'supchan')
      tmpcfg.channel = [tmpcfg.channel cfg.supchan(:)'];
    end

    % select the data in the channels of interest
    [dum, datchanindx] = match_str(tmpcfg.channel, data.label);
    if Ntrials==1
      data.avg = data.avg(datchanindx,:);
      data.cov = data.cov(datchanindx,datchanindx);
    else
      data.avg   = data.avg(datchanindx,:);
      data.cov   = data.cov(:,datchanindx,datchanindx);
      data.trial = data.trial(:,datchanindx,:);
    end
    data.label = data.label(datchanindx);

    if isfield(cfg, 'refchan') && ~isempty(cfg.refchan)
      [dum, refchanindx] = match_str(cfg.refchan, data.label);
    else
      refchanindx = [];
    end
    if isfield(cfg, 'supchan') && ~isempty(cfg.supchan)
      [dum, supchanindx] = match_str(cfg.supchan, data.label);
    else
      supchanindx = [];
    end
    Nchans = length(tmpcfg.channel); % update the number of channels

  else
    % HACK: use the default code
    % select the channels of interest
    [dum, datchanindx] = match_str(cfg.channel, data.label);
    if strcmp(data.dimord, 'chan_time')
      % It is in principle possible to have timelockanalysis with
      % keeptrial=yes and only a single trial in the raw data.
      % In that case the covariance should be represented as Nchan*Nchan
      data.avg = data.avg(datchanindx,:);
      data.cov = reshape(data.cov, length(datchanindx), length(datchanindx));
      data.cov = data.cov(datchanindx,datchanindx);
    else
      data.avg   = data.avg(datchanindx,:);
      data.cov   = data.cov(:,datchanindx,datchanindx);
      data.trial = data.trial(:,datchanindx,:);
    end
    data.label = data.label(datchanindx);
    Nchans     = length(data.label);
  end

  if nargin>2
    % baseline and active are only available together for resampling purposes,
    % hence I assume here that there are multiple trials in both
    baseline.avg   = baseline.avg(datchanindx,:);
    baseline.cov   = baseline.cov(:,datchanindx,datchanindx);
    baseline.trial = baseline.trial(:,datchanindx,:);
    % this is required for averaging 2 conditions using prepare_resampled_data
    cfg2 = [];
    cfg2.numcondition = 2;
  end

  % prepare the resampling of the trials, or average the data if multiple trials are present and no resampling is neccessary
  if (strcmp(cfg.jackknife, 'yes') || strcmp(cfg.bootstrap, 'yes') || strcmp(cfg.pseudovalue, 'yes') || strcmp(cfg.singletrial, 'yes') || strcmp(cfg.rawtrial, 'yes') || strcmp(cfg.randomization, 'yes')) && ~strcmp(data.dimord, 'rpt_chan_time')
    error('multiple trials required in the data\n');

  elseif strcmp(cfg.permutation, 'yes')
    % compute the average and covariance without resampling
    [dum, avgA, covA, avgB, covB] = prepare_resampled_data(cfg2 , data.trial, data.cov, baseline.trial, baseline.cov);
    % compute the average and covariance with random permutation
    [cfg, avRA, coRA, avRB, coRB] = prepare_resampled_data(cfg, data.trial, data.cov, baseline.trial, baseline.cov);
    % concatenate the different resamplings
    avg = cat(1, reshape(avgA, [1 Nchans Nsamples]), reshape(avgB, [1 Nchans Nsamples]), avRA, avRB);
    Cy  = cat(1, reshape(covA, [1 Nchans Nchans  ]), reshape(covB, [1 Nchans Nchans  ]), coRA, coRB);
    % clear temporary working copies
    clear avgA avgB covA covB
    clear avRA avRB coRA coRB
    % the order of the resamplings should be [avgA avgB randA randB randA randB randA randB ....]
    Nrepetitions = 2*cfg.numpermutation + 2;
    order = [1 2 3:2:Nrepetitions 4:2:Nrepetitions];
    avg = avg(order,:,:);
    Cy  = Cy (order,:,:);

  elseif strcmp(cfg.randomization, 'yes')
    % compute the average and covariance without resampling
    [dum, avgA, covA, avgB, covB] = prepare_resampled_data(cfg2 , data.trial, data.cov, baseline.trial, baseline.cov);
    % compute the average and covariance with random resampling
    [cfg, avRA, coRA, avRB, coRB] = prepare_resampled_data(cfg, data.trial, data.cov, baseline.trial, baseline.cov);
    % concatenate the different resamplings
    avg = cat(1, reshape(avgA, [1 Nchans Nsamples]), reshape(avgB, [1 Nchans Nsamples]), avRA, avRB);
    Cy  = cat(1, reshape(covA, [1 Nchans Nchans  ]), reshape(covB, [1 Nchans Nchans  ]), coRA, coRB);
    % clear temporary working copies
    clear avgA avgB covA covB
    clear avRA avRB coRA coRB
    % the order of the resamplings should be [avgA avgB randA randB randA randB randA randB ....]
    Nrepetitions = 2*cfg.numrandomization + 2;
    order = [1 2 3:2:Nrepetitions 4:2:Nrepetitions];
    avg = avg(order,:,:);
    Cy  = Cy (order,:,:);

  elseif strcmp(cfg.jackknife, 'yes')
    % compute the jackknife repetitions for the average and covariance
    [cfg, avg, Cy] = prepare_resampled_data(cfg, data.trial, data.cov);
    Nrepetitions = Ntrials;

  elseif strcmp(cfg.bootstrap, 'yes')
    % compute the bootstrap repetitions for the average and covariance
    [cfg, avg, Cy] = prepare_resampled_data(cfg, data.trial, data.cov);
    Nrepetitions = cfg.numbootstrap;

  elseif strcmp(cfg.pseudovalue, 'yes')
    % compute the pseudovalue repetitions for the average and covariance
    [cfg, avg, Cy] = prepare_resampled_data(cfg, data.trial, data.cov);
    Nrepetitions = Ntrials+1;

  elseif strcmp(cfg.singletrial, 'yes')
    % The idea is that beamformer uses the average covariance to construct the
    % filter and applies it to the single trial covariance/csd. The problem
    % is that beamformer will use the averaged covariance/csd to estimate the
    % power and not the single trial covariance/csd
    error('this option contains a bug, and is therefore not supported at the moment');
    % average the single-trial covariance matrices
    Cy = mean(data.cov,1);
    % copy the average covariance matrix for every individual trial
    Cy = repmat(Cy, [Ntrials 1 1]);
    % keep the single-trial ERFs, rename them to avg for convenience
    avg = data.trial;
    Nrepetitions = Ntrials;

  elseif strcmp(cfg.rawtrial, 'yes')
    % do not do any resampling, keep the single-trial covariances
    Cy = data.cov;
    % do not do any resampling, keep the single-trial ERFs (rename them to avg for convenience)
    avg = data.trial;
    Nrepetitions = Ntrials;

  elseif Ntrials>1
    % average the single-trial covariance matrices
    Cy  = reshape(mean(data.cov,1), [Nchans Nchans]);
    % select the average ERF
    avg = data.avg;
    Nrepetitions = 1;

  elseif Ntrials==1
    % select the average covariance matrix
    Cy  = data.cov;
    % select the average ERF
    avg = data.avg;
    Nrepetitions = 1;
  end

  % reshape so that it also looks like one trial (out of many)
  if Nrepetitions==1
    Cy  = reshape(Cy , [1 Nchans Nchans]);
    avg = reshape(avg, [1 Nchans Nsamples]);
  end

  % get the relevant low level options from the cfg and convert into key-value pairs
  optarg = cfg2keyval(getfield(cfg, cfg.method));

  if strcmp(cfg.method, 'lcmv') && ~isfield(grid, 'filter'),
    for i=1:Nrepetitions
      fprintf('scanning repetition %d\n', i);
      dip(i) = beamformer_lcmv(grid, sens, vol, squeeze(avg(i,:,:)), squeeze(Cy(i,:,:)), optarg{:});
    end
  elseif strcmp(cfg.method, 'lcmv')
    %don't loop over repetitions (slow), but reshape the input data to obtain single trial timecourses efficiently
    %in the presence of filters pre-computed on the average (or whatever)
    siz    = size(avg);
    tmpavg = reshape(permute(avg,[2 3 1]),[siz(2) siz(3)*siz(1)]);
    tmpdip = beamformer_lcmv(grid, sens, vol, tmpavg, squeeze(mean(Cy,1)), optarg{:});
    for i=1:length(tmpdip.inside)
      indx = tmpdip.inside(i);
      tmpdip.mom{indx} = reshape(tmpdip.mom{indx}, [siz(3) siz(1)])';
    end
    try, tmpdip = rmfield(tmpdip, 'pow'); end
    try, tmpdip = rmfield(tmpdip, 'cov'); end
    try, tmpdip = rmfield(tmpdip, 'noise'); end
    for i=1:Nrepetitions
      dip(i).pos     = tmpdip.pos;
      dip(i).inside  = tmpdip.inside;
      dip(i).outside = tmpdip.outside;
      dip(i).mom     = cell(1,size(tmpdip.pos,1));
      for ii=1:length(tmpdip.inside)
        indx = tmpdip.inside(ii);
	dip(i).mom{indx} = tmpdip.mom{indx}(i,:);
      end
    end
  elseif strcmp(cfg.method, 'sam')
    for i=1:Nrepetitions
      fprintf('scanning repetition %d\n', i);
      dip(i) = beamformer_sam(grid, sens, vol, squeeze(avg(i,:,:)), squeeze(Cy(i,:,:)), optarg{:});
    end
  elseif strcmp(cfg.method, 'pcc')
    for i=1:Nrepetitions
      fprintf('scanning repetition %d\n', i);
      dip(i) = beamformer_pcc(grid, sens, vol, squeeze(avg(i,:,:)), squeeze(Cy(i,:,:)), optarg{:}, 'refdip', cfg.refdip, 'refchan', refchanindx, 'supchan', supchanindx);
    end
  elseif strcmp(cfg.method, 'mne')
    for i=1:Nrepetitions
      fprintf('estimating current density distribution for repetition %d\n', i);
      dip(i) = minimumnormestimate(grid, sens, vol, squeeze(avg(i,:,:)),                     optarg{:});
    end
  elseif strcmp(cfg.method, 'loreta')
    for i=1:Nrepetitions
      fprintf('estimating LORETA current density distribution for repetition %d\n', i);
      dip(i) = loreta(             grid, sens, vol, squeeze(avg(i,:,:)),                     optarg{:});
    end
  elseif strcmp(cfg.method, 'rv')
    for i=1:Nrepetitions
      fprintf('estimating residual variance at each grid point for repetition %d\n', i);
      dip(i) = residualvariance(   grid, sens, vol, squeeze(avg(i,:,:)),                     optarg{:});
    end
  elseif strcmp(cfg.method, 'music')
    for i=1:Nrepetitions
      fprintf('computing multiple signal classification for repetition %d\n', i);
      if hascovariance
        dip(i) = music(grid, sens, vol, squeeze(avg(i,:,:)), 'cov', squeeze(Cy(i,:,:)), optarg{:});
      else
        dip(i) = music(grid, sens, vol, squeeze(avg(i,:,:)),                            optarg{:});
      end
    end
  elseif strcmp(cfg.method, 'mvl')
    for i=1:Nrepetitions
      fprintf('estimating current density distribution for repetition %d\n', i);
      fns = fieldnames(cfg);
      optarg = cell(1,length(fns));
      n=1;
      for c=1:length(fns)
        optarg{n} = fns{c};
        optarg{n+1} = cfg.(fns{c});
        n=n+2;
      end
      dip(i) = mvlestimate(grid, sens, vol, squeeze(avg(i,:,:)), optarg{:});
    end
  else
    error(sprintf('method ''%s'' is unsupported for source reconstruction in the time domain', cfg.method));
  end

end % if freq or timelock data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean up and collect the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(grid, 'xgrid')
  % the dipoles are placed on a regular grid that is aligned with the carthesian axes
  % copy the description of the grid axes
  source.xgrid = grid.xgrid;
  source.ygrid = grid.ygrid;
  source.zgrid = grid.zgrid;
  source.dim   = [length(source.xgrid) length(source.ygrid) length(source.zgrid)];
else
  source.dim   = [size(grid.pos,1) 1];
end

source.vol = vol;
if exist('grad', 'var')
  source.grad = grad;
elseif exist('elec', 'var')
  source.elec = elec;
end

if istimelock
  % add the time axis to the output
  source.time = data.time;
elseif iscomp
  % FIXME, add the component numbers to the output
elseif isfreq
  % add the frequency axis to the output
  cfg.frequency    = data.freq(nearest(data.freq, cfg.frequency));
  source.frequency = cfg.frequency;
  if isfield(data, 'time') && isfield(cfg, 'latency')
    cfg.latency    = data.time(nearest(data.time, cfg.latency));
    source.latency = cfg.latency;
  end
  if isfield(data, 'cumtapcnt'),
    source.cumtapcnt = data.cumtapcnt;
  end
end

if exist('dip', 'var')
  % do some cleaning up, keep the dipole positions etc. in the global structure and not in each trial
  source.pos     = dip(1).pos;
  source.inside  = dip(1).inside;
  source.outside = dip(1).outside;
  dip = rmfield(dip, 'pos');
  dip = rmfield(dip, 'inside');
  dip = rmfield(dip, 'outside');
  if isfield(dip(1), 'leadfield')
    source.leadfield = dip(1).leadfield;
    dip = rmfield(dip, 'leadfield');
  end
elseif exist('grid', 'var')
  % no scanning has been done, probably only the leadfield has been computed
  try, source.pos       = grid.pos;       end
  try, source.inside    = grid.inside;    end
  try, source.outside   = grid.outside;   end
  try, source.leadfield = grid.leadfield; end
end

if strcmp(cfg.keepleadfield, 'yes') && ~isfield(source, 'leadfield')
  % add the precomputed leadfields to the output source
  source.leadfield = grid.leadfield;
elseif strcmp(cfg.keepleadfield, 'no') && isfield(source, 'leadfield')
  % remove the precomputed leadfields from the output source
  source = rmfield(source, 'leadfield');
end

% remove the precomputed leadfields from the cfg regardless of what keepleadfield is saying
% it should not be kept in cfg, since there it takes up too much space
if isfield(cfg, 'grid') && isfield(cfg.grid, 'leadfield')
  cfg.grid = rmfield(cfg.grid, 'leadfield');
end

if convertfreq
  % FIXME, convert the source reconstruction back to a frequency representation
elseif convertcomp
  % FIXME, convert the source reconstruction back to a component representation
end

if strcmp(cfg.jackknife, 'yes')
  source.method = 'jackknife';
  source.trial = dip;
  source.df = Ntrials;
elseif strcmp(cfg.bootstrap, 'yes')
  source.method = 'bootstrap';
  source.trial = dip;
  source.df = Ntrials;
elseif strcmp(cfg.pseudovalue, 'yes')
  source.method = 'pseudovalue';
  source.trial = dip;
elseif strcmp(cfg.singletrial, 'yes')
  source.method = 'singletrial';
  source.trial = dip;
  source.df = Ntrials;     % is this correct?
elseif strcmp(cfg.rawtrial, 'yes')
  source.method = 'rawtrial';
  source.trial = dip;
  source.df = Ntrials;     % is this correct?
elseif strcmp(cfg.randomization, 'yes')
  % assign the randomized resamplings to the output, keeping the special meaning of trial 1 and 2 in mind
  source.method = 'randomization';
  source.avgA = dip(1);
  source.avgB = dip(2);
  source.trialA = dip(1+2*(1:cfg.numrandomization));
  source.trialB = dip(2+2*(1:cfg.numrandomization));
elseif strcmp(cfg.permutation, 'yes')
  % assign the randomized resamplings to the output, keeping the special meaning of trial 1 and 2 in mind
  source.method = 'permutation';
  source.avgA = dip(1);
  source.avgB = dip(2);
  source.trialA = dip(1+2*(1:cfg.numpermutation));
  source.trialB = dip(2+2*(1:cfg.numpermutation));
elseif exist('dip', 'var')
  % it looks like beamformer analysis was done on an average input, keep the average source reconstruction
  source.method = 'average';
  source.avg = dip;
else
  % apparently no computations were performed
end

if (strcmp(cfg.jackknife, 'yes') || strcmp(cfg.bootstrap, 'yes') || strcmp(cfg.pseudovalue, 'yes') || strcmp(cfg.singletrial, 'yes') || strcmp(cfg.rawtrial, 'yes')) && strcmp(cfg.keeptrials, 'yes')
  % keep the source reconstruction for each repeated or resampled trial
  source.trial	= dip;
end

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: sourceanalysis.m,v 1.141 2009/10/12 14:44:07 jansch Exp $';
% remember the configuration details of the input data
if nargin==2
  try, cfg.previous    = data.cfg;     end
elseif nargin==3
  cfg.previous = [];
  try, cfg.previous{1} = data.cfg;     end
  try, cfg.previous{2} = baseline.cfg; end
end
% remember the exact configuration details in the output
source.cfg = cfg;

fprintf('total time in sourceanalysis %.1f seconds\n', toc);


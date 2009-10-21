function [scd] = scalpcurrentdensity(cfg, data);

% SCALPCURRENTDENSITY computes an estimate of the SCD using the
% second-order derivative (the surface Laplacian) of the EEG potential
% distribution
%
% Use as
%   [data] = scalpcurrentdensity(cfg, data)
% or
%   [timelock] = scalpcurrentdensity(cfg, timelock)
% where the input data is obtained from PREPROCESSING or from
% TIMELOCKANALYSIS. The output data has the same format as the input
% and can be used in combination with most other FieldTrip functions
% (e.g. FREQNALYSIS or TOPOPLOTER).
%
% The configuration can contain
%   cfg.method       = 'finite' for finite-difference method or
%                      'spline' for spherical spline method
%                      'hjorth' for Hjorth approximation method
%   cfg.elecfile     = string, file containing the electrode definition
%   cfg.elec         = structure with electrode definition
%   cfg.conductivity = conductivity of the skin (default = 0.33 S/m)
%   cfg.trials       = 'all' or a selection given as a 1xN vector (default = 'all')
%
% Note that the skin conductivity, electrode dimensions and the potential
% all have to be expressed in the same SI units, otherwise the units of
% the SCD values are not scaled correctly. The spatial distribution still
% will be correct.
%
% The 'finite' method implements
%   TF Oostendorp, A van Oosterom; The surface Laplacian of the potential:
%   theory and application. IEEE Trans Biomed Eng, 43(4): 394-405, 1996.
%   G Huiskamp; Difference formulas for the surface Laplacian on a
%   triangulated sphere. Journal of Computational Physics, 2(95): 477-496,
%   1991.
%
% The 'spline' method implements
%   F. Perrin, J. Pernier, O. Bertrand, and J. F. Echallier.
%   Spherical splines for scalp potential and curernt density mapping.
%   Electroencephalogr Clin Neurophysiol, 72:184-187, 1989
% including their corrections in
%   F. Perrin, J. Pernier, O. Bertrand, and J. F. Echallier.
%   Corrigenda: EEG 02274, Electroencephalography and Clinical
%   Neurophysiology 76:565.
%
% The 'hjorth' method implements
%   B. Hjort; An on-line transformation of EEG ccalp potentials into
%   orthogonal source derivation. Electroencephalography and Clinical
%   Neurophysiology 39:526-530, 1975.

% Copyright (C) 2004-2006, Robert Oostenveld
%
% $Log: scalpcurrentdensity.m,v $
% Revision 1.25  2008/11/12 19:26:17  roboos
% documented hjorth, added support for layout in constructing of hjorth
%
% Revision 1.24  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.23  2008/07/21 13:19:05  roboos
% implemented Hjorth filtering, using neighbourselection and apply_warp
% changed finite method to use apply_warp
%
% Revision 1.22  2008/05/06 14:23:32  sashae
% change in trial selection, cfg.trials can be a logical
%
% Revision 1.21  2008/04/29 14:35:20  roboos
% fixed bug in looping over trials, thanks to Manuel
%
% Revision 1.20  2008/03/05 10:46:36  roboos
% moved electrode reading functionality from read_fcdc_elec to read_sens, switched to the use of the new function
%
% Revision 1.19  2008/01/31 17:20:17  sashae
% added option for trial selection
%
% Revision 1.18  2007/05/30 07:19:27  roboos
% ensure that the input is not MEG data
%
% Revision 1.17  2007/05/02 15:59:13  roboos
% be more strict on the input and output data: It is now the task of
% the private/checkdata function to convert the input data to raw
% data (i.e. as if it were coming straight from preprocessing).
% Furthermore, the output data is NOT converted back any more to the
% input data, i.e. the output data is the same as what it would be
% on raw data as input, regardless of the actual input.
%
% Revision 1.16  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.15  2007/03/27 11:05:19  ingnie
% changed call to fixdimord in call to checkinput
%
% Revision 1.14  2006/09/07 12:31:00  roboos
% updated documentation
%
% Revision 1.13  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.12  2006/02/23 10:28:16  roboos
% changed dimord strings for consistency, changed toi and foi into time and freq, added fixdimord where neccessary
%
% Revision 1.11  2005/08/05 08:51:14  roboos
% switched to use of read_fcdc_elec subfunction for reading electrodes
% switched to use of data2raw and raw2data subfunctions for converting between avg and raw
%
% Revision 1.10  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.9  2005/06/28 19:50:38  roboos
% minor change to solve problem when compiling to c-code, functionality is the same
% load(cfg.filename) gives compilation error, so first copy the string with the filename from the structure in a plain variable and then load(...)
%
% Revision 1.8  2004/05/05 13:55:20  roberto
% fixed bug that occured when average data was made with keeptrials=yes
%
% Revision 1.7  2004/04/13 16:31:09  roberto
% fixed bug in dbstack selection of function filename for Matlab 6.1
%
% Revision 1.6  2004/04/13 14:25:24  roberto
% wrapped code-snippet around mfilename to make it compatible with Matlab 6.1
%
% Revision 1.5  2004/03/29 15:10:18  roberto
% monor update to handling of average data
%
% Revision 1.4  2004/03/22 15:56:35  roberto
% restructured version information in output configuration
%
% Revision 1.3  2004/03/22 15:43:17  roberto
% attempt to fix CVS version
%
% Revision 1.2  2004/03/22 15:36:00  roberto
% added version information to cfg
%
% Revision 1.1  2004/03/22 15:35:08  roberto
% first version
%

fieldtripdefs

% check if the input data is valid for this function
data = checkdata(data, 'datatype', 'raw', 'feedback', 'yes', 'ismeg', 'no');

% set the defaults
if ~isfield(cfg, 'method'),        cfg.method = 'spline';    end
if ~isfield(cfg, 'conductivity'),  cfg.conductivity = 0.33;  end    % in S/m
if ~isfield(cfg, 'trials'),        cfg.trials = 'all';       end

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  if islogical(cfg.trials),  cfg.trials=find(cfg.trials);  end
  fprintf('selecting %d trials\n', length(cfg.trials));
  data.trial  = data.trial(cfg.trials);
  data.time   = data.time(cfg.trials);
  % update the trial definition (trl)
  if isfield(data, 'cfg') % try to locate the trl in the nested configuration
    trl = findcfg(data.cfg, 'trl');
  else
    trl = [];
  end
  if isempty(trl)
    % a trial definition is expected in each continuous data set
    warning('could not locate the trial definition ''trl'' in the data structure');
  else
    cfg.trlold=trl;
    cfg.trl=trl(cfg.trials,:);
  end
end

% get the electrode positions
if isfield(cfg, 'elecfile')
  fprintf('reading electrodes from file %s\n', cfg.elecfile);
  elec = read_sens(cfg.elecfile);
elseif isfield(cfg, 'elec')
  fprintf('using electrodes specified in the configuration\n');
  elec = cfg.elec;
elseif isfield(data, 'elec')
  fprintf('using electrodes specified in the data\n');
  elec = data.elec;
elseif isfield(cfg, 'layout')
  fprintf('using the 2-D layout to determine the neighbours\n');
  cfg.layout = prepare_layout(cfg);
  cfg.neighbours = neighbourselection(cfg, data);
  % create a dummy electrode structure, this is needed for channel selection
  elec = [];
  elec.label  = cfg.layout.label;
  elec.pnt    = cfg.layout.pos;
  elec.pnt(:,3) = 0;
else
  error('electrode positions were not specified');
end

% remove all junk fields from the electrode array
tmp = elec;
elec = [];
elec.pnt = tmp.pnt;
elec.label = tmp.label;

% find matching electrode positions and channels in the data
[dataindx, elecindx] = match_str(data.label, elec.label);
data.label = data.label(dataindx);
elec.label = elec.label(elecindx);
elec.pnt   = elec.pnt(elecindx, :);
Ntrials = length(data.trial);
for trlop=1:Ntrials
  data.trial{trlop} = data.trial{trlop}(dataindx,:);
end

% compute SCD for each trial
if strcmp(cfg.method, 'spline')
  for trlop=1:Ntrials
    % do not compute intrepolation, but only one value at [0 0 1]
    % this also gives L1, the laplacian of the original data in which we
    % are interested here
    fprintf('computing SCD for trial %d\n', trlop);
    [V2, L2, L1] = splint(elec.pnt, data.trial{trlop}, [0 0 1]);
    scd.trial{trlop} = L1;
  end

elseif strcmp(cfg.method, 'finite')
  % the finite difference approach requires a triangulation
  prj = elproj(elec.pnt);
  tri = delaunay(prj(:,1), prj(:,2));
  % the new electrode montage only needs to be computed once for all trials
  montage.tra = lapcal(elec.pnt, tri);
  montage.labelorg = data.label;
  montage.labelnew = data.label;
  % apply the montage to the data, also update the electrode definition
  scd = apply_montage(data, montage);
  elec = apply_montage(elec, montage);

elseif strcmp(cfg.method, 'hjorth')
  % the Hjorth filter requires a specification of the neighbours
  if ~isfield(cfg, 'neighbours')
    tmpcfg      = [];
    tmpcfg.elec = elec;
    cfg.neighbours = neighbourselection(tmpcfg, data);
  end
  % convert the neighbourhood structure into a montage
  labelnew = {};
  labelorg = {};
  for i=1:length(cfg.neighbours)
    labelnew  = cat(2, labelnew, cfg.neighbours{i}.label);
    labelorg = cat(2, labelorg, cfg.neighbours{i}.neighblabel(:)');
  end
  labelorg = cat(2, labelnew, labelorg);
  labelorg = unique(labelorg);
  tra = zeros(length(labelnew), length(labelorg));
  for i=1:length(cfg.neighbours)
    thischan   = match_str(labelorg, cfg.neighbours{i}.label);
    thisneighb = match_str(labelorg, cfg.neighbours{i}.neighblabel);
    tra(i, thischan) = 1;
    tra(i, thisneighb) = -1/length(thisneighb);
  end
  % combine it in a montage 
  montage.tra = tra;
  montage.labelorg = labelorg;
  montage.labelnew = labelnew;
  % apply the montage to the data, also update the electrode definition
  scd = apply_montage(data, montage);
  elec = apply_montage(elec, montage);

else
  error('unknown method for SCD computation');
end

if strcmp(cfg.method, 'spline') || strcmp(cfg.method, 'finite')
  % correct the units
  warning('trying to correct the units, assuming uV and mm');
  for trlop=1:Ntrials
    % The surface laplacian is proportional to potential divided by squared distance which means that, if
    % - input potential is in uV, which is 10^6 too large
    % - units of electrode positions are in mm, which is 10^3 too large
    % these two cancel out against each other. Hence the computed laplacian
    % is in SI units (MKS).
    scd.trial{trlop} = cfg.conductivity * -1 * scd.trial{trlop};
  end
  fprintf('output surface laplacian is in V/m^2');
else
  fprintf('output Hjorth filtered potential is in uV');
end

% collect the results
scd.elec    = elec;
scd.time    = data.time;
scd.label   = data.label;
scd.fsample = data.fsample;

% store the configuration of this function call, including that of the previous function call
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id   = '$Id: scalpcurrentdensity.m,v 1.25 2008/11/12 19:26:17 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
scd.cfg = cfg;


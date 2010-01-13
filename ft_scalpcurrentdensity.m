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
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

% check if the input data is valid for this function
data = checkdata(data, 'datatype', 'raw', 'feedback', 'yes', 'ismeg', 'no');

% set the defaults
if ~isfield(cfg, 'method'),        cfg.method = 'spline';    end
if ~isfield(cfg, 'conductivity'),  cfg.conductivity = 0.33;  end    % in S/m
if ~isfield(cfg, 'trials'),        cfg.trials = 'all';       end

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  fprintf('selecting %d trials\n', length(cfg.trials));
  data = selectdata(data, 'rpt', cfg.trials);  
  if isfield(data, 'cfg') % try to locate the trl in the nested configuration
    cfg.trlold = findcfg(data.cfg, 'trlold');
    cfg.trl    = findcfg(data.cfg, 'trl');
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
cfg.version.id   = '$Id$';
% remember the configuration details of the input data
try, cfg.previous = data.cfg; end
% remember the exact configuration details in the output
scd.cfg = cfg;


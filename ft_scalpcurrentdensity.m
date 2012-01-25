function [scd] = ft_scalpcurrentdensity(cfg, data)

% FT_SCALPCURRENTDENSITY computes an estimate of the SCD using the
% second-order derivative (the surface Laplacian) of the EEG potential
% distribution
%
% Use as
%   [data] = ft_scalpcurrentdensity(cfg, data)
% or
%   [timelock] = ft_scalpcurrentdensity(cfg, timelock)
% where the input data is obtained from FT_PREPROCESSING or from
% FT_TIMELOCKANALYSIS. The output data has the same format as the input
% and can be used in combination with most other FieldTrip functions
% such as FT_FREQNALYSIS or FT_TOPOPLOTER.
%
% The configuration can contain
%   cfg.method       = 'finite' for finite-difference method or
%                      'spline' for spherical spline method
%                      'hjorth' for Hjorth approximation method
%   cfg.elecfile     = string, file containing the electrode definition
%   cfg.elec         = structure with electrode definition
%   cfg.trials       = 'all' or a selection given as a 1xN vector (default = 'all')
%
% The spline and finite method require the following
%   cfg.conductivity = conductivity of the skin (default = 0.33 S/m)
% 
% The hjorth method requires the following
%   cfg.neighbours   = neighbourhood structure, see FT_PREPARE_NEIGHBOURS
%
% Note that the skin conductivity, electrode dimensions and the potential
% all have to be expressed in the same SI units, otherwise the units of
% the SCD values are not scaled correctly. The spatial distribution still
% will be correct.
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
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
%   B. Hjort; An on-line transformation of EEG scalp potentials into
%   orthogonal source derivation. Electroencephalography and Clinical
%   Neurophysiology 39:526-530, 1975.

% Copyright (C) 2004-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar data

% set the defaults
if ~isfield(cfg, 'method'),        cfg.method = 'spline';    end
if ~isfield(cfg, 'conductivity'),  cfg.conductivity = 0.33;  end    % in S/m
if ~isfield(cfg, 'trials'),        cfg.trials = 'all';       end

if strcmp(cfg.method, 'hjorth')
  cfg = ft_checkconfig(cfg, 'required', {'neighbours'});
else
  cfg = ft_checkconfig(cfg); % perform a simple consistency check
end

% store original datatype
dtype = ft_datatype(data);

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes', 'ismeg', 'no');

% select trials of interest
if ~strcmp(cfg.trials, 'all')
  fprintf('selecting %d trials\n', length(cfg.trials));
  data = ft_selectdata(data, 'rpt', cfg.trials);
end

% get the electrode positions
elec = ft_fetch_sens(cfg, data);

% remove all junk fields from the electrode array
tmp  = elec;
elec = [];
elec.chanpos = tmp.chanpos;
elec.elecpos = tmp.elecpos;
elec.label   = tmp.label;

% find matching electrode positions and channels in the data
[dataindx, elecindx] = match_str(data.label, elec.label);
data.label   = data.label(dataindx);
elec.label   = elec.label(elecindx);
elec.chanpos = elec.chanpos(elecindx, :);
Ntrials = length(data.trial);
for trlop=1:Ntrials
  data.trial{trlop} = data.trial{trlop}(dataindx,:);
end

% compute SCD for each trial
if strcmp(cfg.method, 'spline')
  for trlop=1:Ntrials
    % do not compute interpolation, but only one value at [0 0 1]
    % this also gives L1, the laplacian of the original data in which we
    % are interested here
    fprintf('computing SCD for trial %d\n', trlop);
    [V2, L2, L1] = splint(elec.chanpos, data.trial{trlop}, [0 0 1]);
    scd.trial{trlop} = L1;
  end
  
elseif strcmp(cfg.method, 'finite')
  % the finite difference approach requires a triangulation
  prj = elproj(elec.chanpos);
  tri = delaunay(prj(:,1), prj(:,2));
  % the new electrode montage only needs to be computed once for all trials
  montage.tra = lapcal(elec.chanpos, tri);
  montage.labelorg = data.label;
  montage.labelnew = data.label;
  % apply the montage to the data, also update the electrode definition
  scd  = ft_apply_montage(data, montage);
  elec = ft_apply_montage(elec, montage);
  
elseif strcmp(cfg.method, 'hjorth')
  % convert the neighbourhood structure into a montage
  labelnew = {};
  labelorg = {};
  for i=1:length(cfg.neighbours)
    labelnew  = cat(2, labelnew, cfg.neighbours(i).label);
    labelorg = cat(2, labelorg, cfg.neighbours(i).neighblabel(:)');
  end
  labelorg = cat(2, labelnew, labelorg);
  labelorg = unique(labelorg);
  tra = zeros(length(labelnew), length(labelorg));
  for i=1:length(cfg.neighbours)
    thischan   = match_str(labelorg, cfg.neighbours(i).label);
    thisneighb = match_str(labelorg, cfg.neighbours(i).neighblabel);
    tra(i, thischan) = 1;
    tra(i, thisneighb) = -1/length(thisneighb);
  end
  % combine it in a montage
  montage.tra = tra;
  montage.labelorg = labelorg;
  montage.labelnew = labelnew;
  % apply the montage to the data, also update the electrode definition
  scd  = ft_apply_montage(data, montage);
  elec = ft_apply_montage(elec, montage);
  
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
scd.fsample = 1/(data.time{1}(2) - data.time{1}(1));
if isfield(data, 'sampleinfo')
  scd.sampleinfo = data.sampleinfo;
end
if isfield(data, 'trialinfo')
  scd.trialinfo = data.trialinfo;
end

% convert back to input type if necessary
switch dtype
  case 'timelock'
    scd = ft_checkdata(scd, 'datatype', 'timelock');
  otherwise
    % keep the output as it is
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous data

% rename the output variable to accomodate the savevar postamble
data = scd;

ft_postamble history data
ft_postamble savevar data

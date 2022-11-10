function [scd] = ft_scalpcurrentdensity(cfg, data)

% FT_SCALPCURRENTDENSITY computes an estimate of the SCD using the
% second-order derivative (the surface Laplacian) of the EEG potential
% distribution
%
% The relation between the surface Laplacian and the SCD is explained
% in more detail on http://tinyurl.com/ptovowl.
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
% The configuration should contain
%   cfg.method       = 'finite' for finite-difference method or
%                      'spline' for spherical spline method
%                      'hjorth' for Hjorth approximation method
%   cfg.elec         = structure with electrode positions or filename, see FT_READ_SENS
%   cfg.trials       = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.feedback     = string, 'no', 'text', 'textbar', 'gui' (default = 'text')
%
% The finite method require the following
%   cfg.conductivity = conductivity of the skin (default = 0.33 S/m)
%
% The spline and finite method require the following
%   cfg.conductivity = conductivity of the skin (default = 0.33 S/m)
%   cfg.lambda       = regularization parameter (default = 1e-05)
%   cfg.order        = order of the splines (default = 4)
%   cfg.degree       = degree of legendre polynomials (default for
%                       <=32 electrodes  = 9,
%                       <=64 electrodes  = 14,
%                       <=128 electrodes = 20,
%                       else             = 32
%
% The hjorth method requires the following
%   cfg.neighbours   = neighbourhood structure, see FT_PREPARE_NEIGHBOURS
%
% For the spline method you can specify the following
%   cfg.badchannel      = cell-array, see FT_CHANNELSELECTION for details (default = [])
%
% Note that the skin conductivity, electrode dimensions and the potential
% all have to be expressed in the same SI units, otherwise the units of
% the SCD values are not scaled correctly. The spatial distribution still
% will be correct.
%
% To facilitate data-handling and distributed computing you can use
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
%
% See also FT_PREPROCESSING, FT_TIMELOCKANALYSIS, FT_FREQNALYSIS, FT_TOPOPLOTER.

% Copyright (C) 2004-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'trial'}); % prevent accidental typos, see issue 1729

% set the defaults
cfg.method          = ft_getopt(cfg, 'method',       'spline');
cfg.conductivity    = ft_getopt(cfg, 'conductivity', 0.33); % in S/m
cfg.trials          = ft_getopt(cfg, 'trials',       'all', 1);
cfg.feedback        = ft_getopt(cfg, 'feedback',     'text');
cfg.badchannel      = ft_getopt(cfg, 'badchannel',     {});

switch cfg.method
  case 'hjorth'
    cfg = ft_checkconfig(cfg, 'required', {'neighbours'});
  case 'spline'
    cfg.lambda  = ft_getopt(cfg, 'lambda', 1e-5);
    cfg.order   = ft_getopt(cfg, 'order', 4);
    cfg.degree  = ft_getopt(cfg, 'degree', []);

    if isempty(cfg.degree) % determines degree of Legendre polynomials bases on number of electrodes
      nchan = numel(data.label);
      if nchan<=32
        cfg.degree = 9;
      elseif nchan<=64
        cfg.degree = 14;
      elseif nchan<=128
        cfg.degree = 20;
      else
        cfg.degree = 32;
      end
    end
  otherwise
    cfg = ft_checkconfig(cfg); % perform a simple consistency check
end

% store original datatype
dtype = ft_datatype(data);

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes', 'ismeg', []);

% get the electrode positions
tmpcfg = cfg;
tmpcfg.senstype = 'EEG';
elec = ft_fetch_sens(tmpcfg, data);

% select channels and trials of interest
tmpcfg = keepfields(cfg, {'trials', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
tmpcfg.channel = elec.label;
data = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

Ntrials = numel(data.trial);

if isempty(cfg.badchannel)
  % check if the first sample of the first trial contains NaNs; if so treat it as a bad channel
  cfg.badchannel = ft_channelselection(find(isnan(data.trial{1}(:,1))), data.label);
end

% match the order of the data channels with the channel positions, order them according to the data
[datindx, elecindx] = match_str(data.label, elec.label);
[goodindx, tmp]     = match_str(data.label, setdiff(data.label, cfg.badchannel, 'stable'));

if ~isempty(cfg.badchannel)
  ft_info('detected channel %s as bad\n', cfg.badchannel{:});
  tmpcfg         = [];
  tmpcfg.channel = data.label(goodindx);
  data           = ft_selectdata(tmpcfg, data);
end

allchanpos  = elec.chanpos(elecindx,:); % the position of all channels, ordered according to the data
goodchanpos = allchanpos(goodindx,:);   % the position of good channels

% compute SCD for each trial
if strcmp(cfg.method, 'spline')
  fprintf('Checking spherical fit... ');
  [c, r] = fitsphere(allchanpos);
  d = allchanpos - repmat(c, size(allchanpos,1), 1);
  d = sqrt(sum(d.^2, 2));
  d = mean(abs(d) / r);
  if abs(d-1) > 0.1
    ft_warning('bad spherical fit (residual: %.2f%%). The interpolation will be inaccurate.', 100*(d-1));
  elseif abs(d-1) < 0.01
    fprintf('perfect spherical fit (residual: %.1f%%)\n', 100*(d-1));
  else
    fprintf('good spherical fit (residual: %.1f%%)\n', 100*(d-1));
  end
  % Builds the spatial filter only once.
  fprintf('Calculating the filter to build the SCD.\n');
  [WVo, WLo] = sphsplint(goodchanpos, allchanpos, cfg.order, cfg.degree, cfg.lambda);
  % Creates a montage to apply the spatial filter.
  montage.tra      = WLo;
  montage.labelold = elec.label(elecindx(goodindx));
  montage.labelnew = elec.label(elecindx);
  % Applies the montage to both the data and electrode definition
  scd  = ft_apply_montage(data, montage);
  elec = ft_apply_montage(elec, montage);

elseif strcmp(cfg.method, 'finite')
  if ~isempty(cfg.badchannel)
    ft_error('the method "%s" does not support the specification of bad channels', cfg.method);
  end
  % the finite difference approach requires a triangulation
  prj = elproj(allchanpos);
  tri = delaunay(prj(:,1), prj(:,2));
  % the new electrode montage only needs to be computed once for all trials
  montage.tra = lapcal(allchanpos, tri);
  montage.labelold = data.label;
  montage.labelnew = data.label;
  % apply the montage to the data, also update the electrode definition
  scd  = ft_apply_montage(data, montage);
  elec = ft_apply_montage(elec, montage);

elseif strcmp(cfg.method, 'hjorth')
  if ~isempty(cfg.badchannel)
    ft_error('the method "%s" does not support the specification of bad channels', cfg.method);
  end
  % convert the neighbourhood structure into a montage
  labelnew = {};
  labelold = {};
  for i=1:length(cfg.neighbours)
    labelnew = cat(2, labelnew, cfg.neighbours(i).label);
    labelold = cat(2, labelold, cfg.neighbours(i).neighblabel(:)');
  end
  labelold = cat(2, labelnew, labelold);
  labelold = unique(labelold);
  tra = zeros(length(labelnew), length(labelold));
  for i=1:length(cfg.neighbours)
    thischan   = match_str(labelold, cfg.neighbours(i).label);
    thisneighb = match_str(labelold, cfg.neighbours(i).neighblabel);
    tra(i, thischan) = 1;
    tra(i, thisneighb) = -1/length(thisneighb);
  end
  % combine it in a montage
  montage.tra = tra;
  montage.labelold = labelold;
  montage.labelnew = labelnew;
  % apply the montage to the data, also update the electrode definition
  scd  = ft_apply_montage(data, montage);
  elec = ft_apply_montage(elec, montage);

else
  ft_error('unknown method "%s"', cfg.method);
end

if strcmp(cfg.method, 'spline') || strcmp(cfg.method, 'finite')
  % correct the units
  ft_warning('trying to correct the units, assuming uV and mm');
  for trlop=1:Ntrials
    % The surface laplacian is proportional to potential divided by squared distance which means that, if
    % - input potential is in uV, which is 10^6 too large
    % - units of electrode positions are in mm, which is 10^3 too large
    % these two cancel out against each other. Hence the computed laplacian
    % is in SI units (MKS).
    scd.trial{trlop} = cfg.conductivity * -1 * scd.trial{trlop};
  end
  fprintf('output surface laplacian is in V/m^2\n');
else
  fprintf('output Hjorth filtered potential is in uV\n');
end

% Adds the electrode definition to the data.
scd.elec = elec;

% convert back to input type if necessary
switch dtype
  case 'timelock'
    scd = ft_checkdata(scd, 'datatype', 'timelock');
  otherwise
    % keep the output as it is
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous data

% rename the output variable to accomodate the savevar postamble
data = scd;

ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data

function [data] = ft_denoise_hfc(cfg,data)

% FT_DENOISE_HFC implements harmonic field correction, which models external
% interference on the recordings as a harmonic magnetic field. It is particulaly
% useful for MEG data with low channel numbers, such as OPM data.
%
% The homogenous field correction method implements Tierney et al. (2021) NIMG,
% https://doi.org/10.1016/j.neuroimage.2021.118484.
%
% The harmonic expansion method implements Tierney et al. (2022) NIMG,
% https://doi.org/10.1016/j.neuroimage.2022.119338.
%
% Use as
%   data = ft_denoise_hfc(cfg,data)
%
% Where cfg is a configuration structure that contains:
%   cfg.channel         = channels for HFC (default = 'all')
% 	cfg.order           = number, spherical harmonic order (default = 1)
%                         order = 1 is a homogenous field
%                         order = 2 includes gradients
%                         order = 3 includes quadratic terms, etc.
%   cfg.trials          = which trials do you want to denoise? (default = 'all')
%   cfg.updatesens      = do you want to update sensor info with projector? (default = 'yes')
%   cfg.feedback        = do you want feedback (default = 'yes')
%   cfg.residualcheck   = do you want to check channel residuals (default = 'yes')
%   cfg.residualthresh  = number in pT, what level of residual signal is fine for quality assurance (default = 50)
%
% See also FT_DENOISE_SYNTHETIC, FT_DENOISE_PCA, FT_DENOISE_DSSP, FT_DENOISE_TSP

% Copyright (C) 2021-22, Tim Tierney, George O'Neill, Robert Seymour, Wellcome Centre for Human Neuroimaging, UCL
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
ft_preamble loadvar    data
ft_preamble provenance data

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% check the input data
data = ft_checkdata(data, 'datatype', {'raw'}, 'ismeg', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'trial'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels'}); % prevent accidental typos, see issue 1729

% set the defaults
cfg.order           = ft_getopt(cfg, 'order', 1);
cfg.channel         = ft_getopt(cfg, 'channel', 'all', 1);
cfg.trials          = ft_getopt(cfg, 'trials', 'all', 1);
cfg.updatesens      = ft_getopt(cfg, 'updatesens', 'yes');
cfg.feedback        = ft_getopt(cfg, 'feedback', 'yes');
cfg.residualcheck   = ft_getopt(cfg, 'residualcheck', 'yes');
cfg.residualthresh  = ft_getopt(cfg, 'residualthresh', 50);

% select trials and channels of interest
tmpcfg = keepfields(cfg, {'trials', 'channel', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo', 'checksize'});
data   = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

% Check match between input data chans and grad chans
[sel, dum] = ismember(data.grad.label, data.label);
num_mismatch = sum(sel == 0);

% Warn the user if there are chans in grad not in data
if num_mismatch > 0
  ft_notice('Found %d channels in grad structure not present in cfg.channel. These will NOT be used for Harmonic Field Correction', num_mismatch);
end

% Check for Dr. Tim Tierney's OPM toolbox on the path, and add if needed
% See: https://github.com/tierneytim/OPM
ft_hastoolbox('opm', 1);

% generate harmonic basis set
opt = [];
opt.li = cfg.order;
opt.v = data.grad.coilpos(sel,:);
opt.o = data.grad.coilori(sel,:);
N = spm_opm_vslm(opt);

% Make montage for the next step
montage             = [];
montage.tra         = eye(length(N)) - N*pinv(N);
montage.labelold    = data.grad.label(sel);
montage.labelnew    = data.grad.label(sel);
montage.chantypeold = data.grad.chantype(sel);
montage.chantypenew = data.grad.chantype(sel);
montage.chanunitold = data.grad.chanunit(sel);
montage.chanunitnew = data.grad.chanunit(sel);

% this is prior to the montage
labelold = data.label;

% apply the montage
data = ft_apply_montage(data, montage, 'keepunused', 'yes');

% Tell the user
ft_info('Applied HFC to the data');

% Update the tra, it is essential to correct the leadfields going forward
if istrue(cfg.updatesens)
  data.grad = ft_apply_montage(data.grad, montage, 'keepunused', 'yes', 'balancename', 'hfc','warning',false);
  ft_info('Converted the sensor description to HFC');
end

% reorder the channels to stay close to the original ordering
[selold, selnew] = match_str(montage.labelold, data.label);
if numel(selnew)==numel(labelold)
  for i=1:numel(data.trial)
    data.trial{i} = data.trial{i}(selnew,:);
  end
  data.label = data.label(selnew);
else
  ft_warning('channel ordering might have changed');
end

% Perform running variance check to identify odd channels
if strcmp(cfg.residualcheck, 'yes')
  residual_check(cfg.residualthresh, data, montage.labelold)
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   data
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function residual_check(residualthresh, data, oldlabels)

% find corrected channels in the output data
[selold2, selnew2] = match_str(oldlabels, data.label);

trvar = [];
for ii = 1:numel(data.trial)
  tmp = data.trial{ii}(selnew2,:);
  Mk  = tmp(:,1);
  Sk  = zeros(size(Mk));
  count = 1;
  for jj = 1:size(tmp,2)
    Xk     = tmp(:,jj);
    Mkprev = Mk;
    Mk     = Mkprev +(Xk-Mkprev)/count;
    Sk     = Sk+(Xk-Mkprev).*(Xk-Mk);
    count=count+1;
  end
  trvar(:,ii) = Sk/(count-1);
end

% Identify the most common chanunit
chanunit = ft_chanunit(data);
[s, i, j] = unique(chanunit);
chanunit = s{mode(j)};

% Express results in pT
scale = ft_scalingfactor(chanunit, 'pT');

SD = mean(sqrt(trvar),2) * scale;

fprintf('Checking for unsual channels post-corrections\n')
count = 0;
for ii = 1:length(SD)
  index = selnew2(ii);
  if SD(ii) > residualthresh
    count = count + 1;
    fprintf('Residual on channel %d (%s) = %3.2f pT\n', index, data.label{index}, SD(ii));
  end
end

if ~count
  fprintf('No unusual channel residuals found!\n')
end

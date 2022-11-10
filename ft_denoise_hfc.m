function [data] = ft_denoise_hfc(cfg,data)

% FT_DENOISE_HFC implements harmonic field correction, which models
% interference as a harmonic magnetic field. It is particulaly useful
% for MEG data with low channel numbers (e.g. OPM data). Homogenous field based
% on Tierney et al. (2021) NIMG, 118484. Harmonic expansion based on Tierney
% et al. (2022) NIMG, 119338.
%
% Use as:
%   data = ft_denoise_hfc(cfg,data)
%
% Where cfg is a configuration structure that contains:
%   cfg.channel         = channels for HFC (default = 'all')
% 	cfg.order           = spherical harmonic order:
%                           order = 1 is a homogenous field; order = 2
%                           includes gradients; order = 3 includes
%                           quadratic terms etc. (default = 1)
%   cfg.trials          = which trials do you want to denoise?
%                           (default = 'all')
%   cfg.updatesens      = do you want to update sensor info with projector?
%                           (default = 'yes')
%   cfg.feedback        = do you want feedback (default = 'yes')
%   cfg.residualcheck   = do you want to check channel residuals
%                          (default = 'yes')
%   cfg.residualthresh  = (in pT) what level of residual signal is fine for
%                           quality assurance (default = 50)
% See also FT_DENOISE_SYNTHETIC, FT_DENOISE_PCA, FT_DENOISE_DSSP, FT_DENOISE_TSP

% Copyright (C) 2021-22, Tim Tierney, George O'Neill, Robert Seymour
% Wellcome Centre for Human Neuroimaging, UCL
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

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'trial'}); % prevent accidental typos, see issue 1729

% check the input data
data = ft_checkdata(data, 'datatype', {'raw'});

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden',  {'trial'}); % prevent accidental typos, see issue 1729
cfg = ft_checkconfig(cfg, 'forbidden',  {'channels'}); % prevent accidental typos, see issue 1729

% Check if the data has a grad structure
if ~isfield(data,'grad')
  error(['Data needs a grad structure']);
end

% set the defaults
cfg.order           = ft_getopt(cfg, 'order', 1);
cfg.channel         = ft_getopt(cfg, 'channel', 'all', 1);
cfg.trials          = ft_getopt(cfg, 'trials', 'all', 1);
cfg.updatesens      = ft_getopt(cfg, 'updatesens', 'yes');
cfg.feedback        = ft_getopt(cfg, 'feedback', 'yes');
cfg.residualcheck   = ft_getopt(cfg, 'residualcheck', 'yes');
cfg.residualthresh  = ft_getopt(cfg, 'residualthresh', 50);

% select data based on cfg.channel and cfg.trials
tmpcfg = keepfields(cfg, {'trials', 'showcallinfo'});

% Select the correct channels using ft_channelselection
tmpcfg.channel = cfg.channel;
data   = ft_selectdata(tmpcfg, data);

% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

% Check match between input data chans and grad chans
[x, dummy] = ismember(data.grad.label,data.label);

% Warn the user if there are chans in grad not in data
num_mismatch = sum(x(:) == 0);

if num_mismatch > 0 && istrue(cfg.feedback)
  ft_warning(['Found ' num2str(num_mismatch) ' channels in grad structure'...
    ' not present in cfg.channel. These channels will NOT be used for '...
    'Harmonic Field Correction']);
end

% Check for Dr. Tim Tierney's OPM toolbox on the path, and add if needed
% See: https://github.com/tierneytim/OPM
ft_hastoolbox('opm', 1);

% generate harmonic basis set
opt = [];
opt.li = cfg.order;
opt.v = data.grad.coilpos(x,:);
opt.o = data.grad.coilori(x,:);
N = spm_opm_vslm(opt);

% Make montage for the next step
montage     = [];
montage.tra = eye(length(N)) - N*pinv(N);
montage.labelold = data.grad.label(x);
montage.labelnew = data.grad.label(x);
montage.chantypeold = data.grad.chantype(x);
montage.chantypenew = data.grad.chantype(x);
montage.chanunitold = data.grad.chanunit(x);
montage.chanunitnew = data.grad.chanunit(x);

labelold = data.label;

% apply montage;
% Tell the user
data = ft_apply_montage(data,montage,'keepunused','yes');
if istrue(cfg.feedback)
  disp('Applied HFC to the data');
end

% Update the tra to account for the g. Essential to correct the lead
% fields going forward.
if istrue(cfg.updatesens)
  data.grad = ft_apply_montage(data.grad, montage, 'keepunused',...
    'yes', 'balancename', 'hfc','warning',false);
  if istrue(cfg.feedback)
    disp('Converted the sensor description to HFC');
  end
end

% reorder the channels to stay close to the original ordering
[dummy, selnew] = match_str(montage.labelold, data.label);
if numel(selnew)==numel(labelold)
  for i=1:numel(data.trial)
    data.trial{i} = data.trial{i}(selnew,:);
  end
  data.label = data.label(selnew);
else
  ft_warning('channel ordering might have changed');
end

% Perform running variance check to identify odd channels
if strcmp(cfg.residualcheck,'yes')
  residual_check(cfg.residualthresh,data,montage.labelold)
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   data
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data


function residual_check(residualthresh,data,oldlabels)
% Script to determine residual variance post HFC

% find corrected channels in the output data
[dummy, selnew2] = match_str(oldlabels, data.label);

trvar = [];
for ii = 1:numel(data.trial)
  tmp = data.trial{ii}(selnew2,:);
  Mk =  tmp(:,1);
  Sk = zeros(size(Mk));
  count = 1;
  for jj = 1:size(tmp,2)
    Xk = tmp(:,jj);
    Mkprev = Mk;
    Mk = Mkprev +(Xk-Mkprev)/count;
    Sk=Sk+(Xk-Mkprev).*(Xk-Mk) ;
    count=count+1;
  end
  trvar(:,ii)=Sk/(count-1);
end

% Identify the most common chanunit
chanunit = ft_chanunit(data);
[s,dummy,j]=unique(chanunit);
chanunit = s{mode(j)};

switch chanunit % some of this are silly, but safety first!
  case 'fT'
    scale = 1e-3;
  case 'pT'
    scale = 1;
  case 'nT'
    scale = 1e3;
  case 'uT'
    scale = 1e6;
  case 'mT'
    scale = 1e9;
  case {'T','T/m'}
    scale = 1e12;
  otherwise
    ft_error('Cannot check residuals due to unknown sensor units!')
end

SD = mean(sqrt(trvar),2)*scale;

fprintf('Checking for unsual channels post-corrections\n')
count = 0;
for ii = 1:length(SD)
  index = selnew2(ii);
  if SD(ii) > residualthresh
    count = count + 1;
    if strcmp(chanunit,'fT')
      fprintf(['Residual on channel ' num2str(index) ', '...
        data.label{index} ': %3.2f pT\n'], SD(ii));
    else
      fprintf(['Residual on channel ' num2str(index) ', '...
        data.label{index} ': %3.2f' chanunit '\n'], SD(ii));
    end
  end
end
if ~count
  fprintf('No unusual channel residuals found!\n')
end

function [data] = ft_denoise_ssp(cfg, varargin)

% FT_DENOISE_SSP projects out topographies based on ambient noise on
% Neuromag/Elekta/MEGIN systems. These topographies are estimated during maintenance
% visits from the engineers of MEGIN. Alternatively, it computes the projectors from
% reference data (e.g., empty room) if it is given as an additional input. For the 
% best results, make sure to preprocess the reference data the same as the data to denoise.
%
% Use as
%   [data] = ft_denoise_ssp(cfg, data)
% or
%   [data] = ft_denoise_ssp(cfg, data, refdata)
% where the input data should come from FT_PREPROCESSING or
% FT_TIMELOCKANALYSIS and the configuration should contain
%   cfg.channel    = the channels to be denoised (default = 'all')
%   cfg.trials     = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.ssp        = 'all' or a cell array of SSP names to apply (default = 'all')
%   cfg.updatesens = 'yes' or 'no', whether to update the sensor array with the spatial projector (default = 'yes')
%
% If refdata is specified, this data will be used for the computation of the spatial projectors, by performing a PCA 
% on-the-fly. In this situation, the configuration should also contain
%   cfg.numcomponent = number of principal components to project out of the data
%                      (default = 3)
%   cfg.refchannel   = the channels used as reference signal (default = 'MEG')
%   cfg.reftrials    = 'all' or a selection given as a 1xN vector, specifies which trials to use 
%                      from the reference data (default = 'all')
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_PREPROCESSING, FT_DENOISE_AMM, FT_DENOISE_DSSP,
% FT_DENOISE_HFC, FT_DENOISE_PCA, FT_DENOISE_PREWHITEN, FT_DENOISE_SSS,
% FT_DENOISE_SYNTHETIC, FT_DENOISE_TSR

% Copyright (C) 2004-2026, Gianpaolo Demarchi, Lau Møller Andersen, Robert Oostenveld, Jan-Mathijs Schoffelen
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
cfg.ssp          = ft_getopt(cfg, 'ssp',          'all');
cfg.trials       = ft_getopt(cfg, 'trials',       'all', 1);
cfg.reftrials    = ft_getopt(cfg, 'reftrials',    'all');
cfg.updatesens   = ft_getopt(cfg, 'updatesens',   'yes');
cfg.numcomponent = ft_getopt(cfg, 'numcomponent', 3);
cfg.component    = ft_getopt(cfg, 'component');
cfg.channel      = ft_getopt(cfg, 'channel',      'all');
cfg.refchannel   = ft_getopt(cfg, 'refchannel',   'MEG');

% check if the input data is valid for this function
if isscalar(varargin)
  data    = ft_checkdata(varargin{1}, 'datatype', 'raw');
  refdata = [];
elseif numel(varargin) == 2
  data    = ft_checkdata(varargin{1}, 'datatype', 'raw', 'feedback', 'yes', 'hassampleinfo', 'yes');
  refdata = ft_checkdata(varargin{2}, 'datatype', 'raw', 'feedback', 'yes', 'hassampleinfo', 'yes');
else
  error('Incorrect number of input arguments.')
end

% store the original type of the input data
dtype = ft_datatype(data);

% check whether it is neuromag data
if isempty(refdata) && ~ft_senstype(data, 'neuromag')
  ft_warning('without additional reference data in the input this function is designed for neuromag data');
end

% select channels and trials of interest
tmpcfg = keepfields(cfg, {'channel', 'trials', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo', 'checksize'});
data   = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

if ~isempty(refdata)
  tmpcfg         = keepfields(cfg, {'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo', 'checksize'});
  tmpcfg.channel = cfg.refchannel;
  tmpcfg.trials  = cfg.reftrials;
  refdata        = ft_selectdata(tmpcfg, refdata);

  % ssp projectors will be computed based on the refdata, requires adjustment of cfg
  cfg.ssp = {'ssp'};
end

% keep track of the original order of the channels
labelold = data.label;

% keep track of the original grad structure
gradorig = data.grad;

if ~isempty(refdata)
  selcomp = 1:cfg.numcomponent;

  ft_info('computing the "ssp" projector\n');
  % compute numcomponent principal components in the reference data
  coeff = pca(cell2mat(refdata.trial)','NumComponents',max(selcomp));
  % compute projector and define montage
  data.grad.balance.ssp.tra = eye(size(coeff,1))-coeff(:,selcomp)*transpose(coeff(:,selcomp));
  data.grad.balance.ssp.labelold = refdata.label;
  data.grad.balance.ssp.labelnew = refdata.label;
end

% first undo/invert the previously applied balancing
while ~isempty(data.grad.balance.current)
  bname    = data.grad.balance.current{end};
  montage = ft_inverse_montage(data.grad.balance.(bname));
  fprintf('reverting the "%s" projection\n', bname);
  data      = ft_apply_montage(data,      montage, 'keepunused', 'yes');
  data.grad = ft_apply_montage(data.grad, montage, 'keepunused', 'no');
  data.grad.balance.current = data.grad.balance.current(1:end-1); % remove this from the list

  if strcmp(bname, 'planar')
    if isfield(data.grad, 'type') && ~isempty(strfind(data.grad.type, '_planar'))
      % remove the _planar postfix from the MEG sensor type
      data.grad.type = sens.type(1:(end-7));
    end
  end
end

% it should be a cell-array
if ischar(cfg.ssp)
  cfg.ssp = {cfg.ssp};
end

if isequal(cfg.ssp, {'all'})
  cfg.ssp = setdiff(fieldnames(data.grad.balance), {'current'});
elseif isequal(cfg.ssp, 'none')
  cfg.ssp = {};
end

% then apply the desired balancing
for i=1:numel(cfg.ssp)
  bname   = cfg.ssp{i};
  montage = data.grad.balance.(bname);
  fprintf('applying the "%s" projection\n', bname);
  data = ft_apply_montage(data, montage, 'keepunused', 'yes');
  if istrue(cfg.updatesens)
    data.grad = ft_apply_montage(data.grad, montage, 'keepunused', 'no');
    data.grad.balance.current{end+1} = bname;
  end
end

% reorder the channels to stay close to the original ordering
[selold, selnew] = match_str(labelold, data.label);
if numel(selnew)==numel(labelold)
  for i=1:numel(data.trial)
    data.trial{i} = data.trial{i}(selnew,:);
  end
  data.label = data.label(selnew);
else
  ft_warning('channel ordering might have changed');
end

% convert back to input type if necessary
switch dtype
  case 'timelock'
    data = ft_checkdata(data, 'datatype', 'timelock');
  otherwise
    % keep the output as it is
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous   data
ft_postamble provenance data
ft_postamble history    data
ft_postamble savevar    data

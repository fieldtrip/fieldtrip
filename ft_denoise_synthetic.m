function [data] = ft_denoise_synthetic(cfg, data)

% FT_DENOISE_SYNTHETIC computes CTF higher-order synthetic gradients for
% preprocessed data and for the corresponding gradiometer definition.
%
% Use as
%   [data] = ft_denoise_synthetic(cfg, data)
% where the input data should come from FT_PREPROCESSING or
% FT_TIMELOCKANALYSIS and the configuration should contain
%   cfg.gradient   = 'none', 'G1BR', 'G2BR' or 'G3BR' specifies the gradiometer
%                    type to which the data should be changed
%   cfg.trials     = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.updatesens = 'no' or 'yes' (default = 'yes')
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
% FT_DENOISE_HFC, FT_DENOISE_PCA, FT_DENOISE_PREWHITEN, FT_DENOISE_SSP,
% FT_DENOISE_SSS, FT_DENOISE_TSR

% Copyright (C) 2004-2025, Robert Oostenveld
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
cfg = ft_checkconfig(cfg, 'required',   {'gradient'});

% set the defaults
cfg.trials     = ft_getopt(cfg, 'trials', 'all', 1);
cfg.updatesens = ft_getopt(cfg, 'updatesens', 'yes');

% store the original type of the input data
dtype = ft_datatype(data);

% check if the input data is valid for this function
% this will convert timelocked input data to a raw data representation if needed
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes', 'hassampleinfo', 'yes');

% check whether it is CTF data
if ~ft_senstype(data, 'ctf')
  ft_error('synthetic gradients can only be computed for CTF data');
end

% check whether there are reference channels in the input data
hasref = ~isempty(ft_channelselection('MEGREF', data.label));
if ~hasref
  ft_error('synthetic gradients can only be computed when the input data contains reference channels');
end

% select trials of interest
tmpcfg = keepfields(cfg, {'trials', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo', 'checksize'});
data   = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

% remember the original channel ordering
labelold = data.label;


% first undo/invert the previously applied balancing
while ~isempty(data.grad.balance.current)
  this_name    = data.grad.balance.current{end};
  this_montage = ft_inverse_montage(data.grad.balance.(this_name));
  fprintf('reverting the "%s" projection\n', this_name);
  data      = ft_apply_montage(data,      this_montage, 'keepunused', 'yes');
  data.grad = ft_apply_montage(data.grad, this_montage, 'keepunused', 'no');
  data.grad.balance.current = data.grad.balance.current(1:end-1); % remove this from the list

  if strcmp(this_name, 'planar')
    if isfield(data.grad, 'type') && ~isempty(strfind(data.grad.type, '_planar'))
      % remove the _planar postfix from the MEG sensor type
      data.grad.type = sens.type(1:(end-7));
    end
  end
end

% then apply the desired balancing
if ~strcmp(cfg.gradient, 'none')
  bname = cfg.gradient;
  montage = data.grad.balance.(bname);
  fprintf('applying the "%s" projection\n', bname);
  data      = ft_apply_montage(data,      montage, 'keepunused', 'yes');
  data.grad = ft_apply_montage(data.grad, montage, 'keepunused', 'no');
  data.grad.balance.current{end+1} = bname;
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

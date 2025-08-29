function [data] = ft_denoise_ssp(cfg, varargin)

% FT_DENOISE_SSP projects out topographies based on ambient noise on
% Neuromag/Elekta/MEGIN systems. These topographies are estimated during maintenance
% visits from the engineers of MEGIN.
% Alternatively, computes projectors from reference data (e.g., empty room) if it 
% is given as an additional input. For best results, make sure to preprocess
% the reference data the same as the data to denoise. 
%
% Use as
%   [data] = ft_denoise_ssp(cfg, data) OR
%   [data] = ft_denoise_ssp(cfg, data, refdata) 
% where data should come from FT_PREPROCESSING and the configuration
% should contain
%   cfg.ssp        = 'all' or a cell array of SSP names to apply (default = 'all')
%   cfg.trials     = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.updatesens = 'no' or 'yes' (default = 'yes')
%   cfg.channel      = Nx1 cell-array with selection of channels from the refdata 
%                      (default = 'all'), see FT_CHANNELSELECTION for details
%
% If refdata is specified the configuration should also contain
%   cfg.numcomponent = number of principal components to project out of the data 
%                      (default = 3)
%
% To facilitate data-handling and distributed cmputing you can use
%   cfg.inputfile   =  ...  
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a sinegle variable, corresponding with the
% input/output structure.
%
% See also FT_PREPROCESSING, FT_DENOISE_SYNTHETIC, FT_DENOISE_PCA

% Copyright (C) 2004-2022, Gianpaolo Demarchi, Lau
% Møller Andersen, Robert Oostenveld, Jan-Mathijs Schoffelen
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
cfg = ft_checkconfig(cfg, 'required',   {'ssp'});

% set the defaults
cfg.ssp        = ft_getopt(cfg, 'ssp', 'all');
cfg.trials     = ft_getopt(cfg, 'trials', 'all', 1);
cfg.updatesens = ft_getopt(cfg, 'updatesens', 'yes');
cfg.numcomponent = ft_getopt(cfg, 'numcomponent', 3);
cfg.channel    = ft_getopt(cfg, 'channel', 'all');

if numel(varargin) == 1
    data = varargin{1};
    refdata = [];
elseif numel(varargin) == 2
    data = varargin{1};
    refdata = varargin{2};
else
    error('Incorrect number of input arguments.')
end

% store the original type of the input data
dtype = ft_datatype(data);

% check if the input data is valid for this function
% this will convert timelocked input data to a raw data representation if needed
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes', 'hassampleinfo', 'yes');

% check whether it is neuromag data
if isempty(refdata) && ~ft_senstype(data, 'neuromag')
  ft_error('SSP without reference data can only be applied to neuromag data');
end

% select trials of interest
tmpcfg = keepfields(cfg, {'channel', 'trials', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo', 'checksize'});
data   = ft_selectdata(tmpcfg, data);

if ~isempty(refdata)
    datachs = ft_channelselection('MEG',data.label); % pick MEG channels from data
    refdata = ft_checkdata(refdata, 'datatype', 'raw', 'feedback', 'yes', 'hassampleinfo', 'yes');
    if ~isempty(setdiff(datachs,refdata.label))
        ft_error('Refdata is missing channels present in data.')
    elseif ~isempty(setdiff(refdata.label,datachs))
        ft_warning('Data is missing channels present in refdata. Dropping extra channels from refdata.')
    end
    cfg.channel = datachs;  
    tmpcfg = keepfields(cfg, {'channel', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo', 'checksize'});
    refdata = ft_selectdata(tmpcfg, refdata);
end

% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

if ~isempty(refdata)
    ft_info('Computing the projector(s)\n');
    % compute numcomponent principal components in the reference data
    [coeff,~,~,~,~] = pca(cell2mat(refdata.trial)','NumComponents',cfg.numcomponent);
    % compute projector and define montage
    data.grad.balance.ssp.tra = eye(size(coeff,1))-coeff*transpose(coeff);
    data.grad.balance.ssp.labelold = refdata.label;
    data.grad.balance.ssp.labelnew = refdata.label;
    if ~isempty(cfg.ssp) && ~isequal(cfg.ssp,{'ssp'})
        ft_error('cfg.ssp was configured incorrectly.');
    else
        cfg.ssp = {'ssp'};
    end
end

% remember the original channel ordering
labelold = data.label;

% apply the balancing to the MEG data and to the gradiometer definition
current = data.grad.balance.current;
desired = cfg.ssp;

if ~strcmp(current, 'none')
  % first undo/invert the previously applied balancing
  try
    current_montage = data.grad.balance.(current);
  catch
    ft_error('unknown balancing for input data');
  end
  fprintf('converting the data from "%s" to "none"\n', current);
  data = ft_apply_montage(data, current_montage, 'keepunused', 'yes', 'inverse', 'yes');
  if istrue(cfg.updatesens)
    fprintf('converting the sensor description from "%s" to "none"\n', current);
    data.grad = ft_apply_montage(data.grad, current_montage, 'keepunused', 'yes', 'inverse', 'yes');
    data.grad.balance.current = 'none';
  end
end % if current

if ~strcmp(desired, 'none')
  % then apply the desired balancing
  if strcmp(desired, 'all')
    desireds = fieldnames(data.grad.balance);
  else
      desireds = cfg.ssp;
    if ~iscell(desireds)
        ft_error('cfg.ssp must be a cell array of projector names')
    end

  end
  for desired_index = 1:length(desireds)
    desired = desireds{desired_index};
    if ~strcmp(desired, 'current')
      try
        desired_montage = data.grad.balance.(desired);
      catch
        ft_error('unknown balancing for input data');
      end
      fprintf('converting the data from "none" to "%s"\n', desired);
      data = ft_apply_montage(data, desired_montage, 'keepunused', 'yes', 'balancename', desired);
      if istrue(cfg.updatesens)
        fprintf('converting the sensor description from "none" to "%s"\n', desired);
        data.grad.balance = rmfield(data.grad.balance,desired); % remove to avoid creating the same balance twice
        data.grad = ft_apply_montage(data.grad, desired_montage, 'keepunused', 'no', 'balancename', desired);
      end
    end % if desired
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

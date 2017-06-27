function [dataout] = ft_channelnormalise(cfg, data)

% FT_CHANNELNORMALISE shifts and scales all channles of the the input data
% to a mean of zero and a standard deviation of one.
%
% Use as
%   [dataout] = ft_channelnormalise(cfg, data)
%
% The configuration can contain
%   cfg.trials = 'all' or a selection given as a 1xN vector (default = 'all')
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_COMPONENTANALYSIS, FT_FREQBASELINE, FT_TIMELOCKBASELINE
%
% Copyright (C) 2010, Jan-Mathijs Schoffelen

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
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set the defaults
cfg.trials = ft_getopt(cfg, 'trials', 'all', 1);

% store original datatype
dtype = ft_datatype(data);

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', 'raw', 'feedback', 'yes');

% select trials of interest
tmpcfg = keepfields(cfg, {'trials', 'showcallinfo'});
data   = ft_selectdata(tmpcfg, data);
% restore the provenance information
[cfg, data] = rollback_provenance(cfg, data);

% initialise some variables
nchan  = numel(data.label);
ntrl   = numel(data.trial);
datsum = zeros(nchan,1);
datssq = zeros(nchan,1);

% create output data, omitting sensor information
% FIXME this can be kept, provided the scaling is built in appropriately
dataout         = [];
dataout.label   = data.label;
dataout.trial   = cell(1,ntrl);
dataout.time    = data.time;

% some fields from the input should be copied over in the output
copyfield = {'fsample', 'sampleinfo', 'trialinfo'};
for i=1:length(copyfield)
  if isfield(data, copyfield{i})
    dataout.(copyfield{i}) = data.(copyfield{i});
  end
end

% compute the mean and std
for k = 1:ntrl
    n(k,1) = size(data.trial{k},2);
    datsum = datsum + sum(data.trial{k},2);
    datssq = datssq + sum(data.trial{k}.^2,2);
end
datmean = datsum./sum(n);
datstd  = sqrt( (datssq - (datsum.^2)./sum(n))./sum(n)); %quick way to compute std from sum and sum-of-squared values

% demean and normalise
for k = 1:ntrl
  dataout.trial{k} = (data.trial{k}-datmean(:,ones(1,n(k))))./datstd(:,ones(1,n(k)));
end

% convert back to input type if necessary
switch dtype
  case 'timelock'
    dataout = ft_checkdata(dataout, 'datatype', 'timelock');
  otherwise
    % keep the output as it is
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble provenance dataout
ft_postamble history    dataout
ft_postamble savevar    dataout

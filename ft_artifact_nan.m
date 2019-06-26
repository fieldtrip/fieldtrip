function [cfg, artifact] = ft_artifact_nan(cfg, data)

% FT_ARTIFACT_NAN identifies artifacts that are indicated in the data as nan (not a
% number) values.
%
% Use as
%   [cfg, artifact] = ft_artifact_nan(cfg, data)
% where the input data is a structure as obtained from FT_REJECTARTIFACT with
% the option cfg.artfctdef.reject='nan'.
%
% The configuration can contain
%   cfg.artfctdef.nan.channel = Nx1 cell-array with selection of channels, see FT_CHANNELSELECTION for details
%
% The output argument "artifact" is a Nx2 matrix comparable to the
% "trl" matrix of FT_DEFINETRIAL. The first column of which specifying the
% beginsamples of an artifact period, the second column contains the
% endsamples of the artifactperiods.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_REJECTARTIFACT, FT_ARTIFACT_CLIP, FT_ARTIFACT_ECG, FT_ARTIFACT_EOG,
% FT_ARTIFACT_JUMP, FT_ARTIFACT_MUSCLE, FT_ARTIFACT_THRESHOLD, FT_ARTIFACT_ZVALUE

% Copyright (C) 2017, Robert Oostenveld
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
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set the default options
cfg.artfctdef             = ft_getopt(cfg, 'artfctdef');
cfg.artfctdef.nan         = ft_getopt(cfg.artfctdef, 'nan');
cfg.artfctdef.nan.channel = ft_getopt(cfg.artfctdef.nan, 'channel', {'all'});

% check if the input data is valid for this function, the input data must be raw
data = ft_checkdata(data, 'datatype', 'raw', 'hassampleinfo', 'yes');

cfg.artfctdef.nan.channel = ft_channelselection(cfg.artfctdef.nan.channel, data.label);
chansel = match_str(data.label, cfg.artfctdef.nan.channel);

artifact = zeros(0,2);

for i=1:numel(data.trial)
  tmp = any(isnan(data.trial{i}(chansel,:)),1);
  if any(tmp)
    % there can be multiple segments with nans
    begsample = find(diff([0 tmp])>0);
    endsample = find(diff([tmp 0])<0);
    for j=1:numel(begsample)
      artifact(end+1,:) = [begsample(j) endsample(j)] + data.sampleinfo(i,1) - 1;
    end
  end
end % for each trial

% remember the details that were used here
cfg.artfctdef.nan          = [];
cfg.artfctdef.nan.artifact = artifact;

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   data
ft_postamble savevar

function [montage, cfg] = ft_prepare_montage(cfg, data)

% FT_PREPARE_MONTAGE creates a referencing scheme based on the input configuration
% options and the channels in the data structure. The resulting montage can be
% given as input to ft_apply_montage, or as cfg.montage to ft_preprocessing.
%
% Use as
%   montage = ft_prepare_montage(cfg, data)
%
% The configuration can contain the following fields:
%   cfg.implicitref   = 'label' or empty, add the implicit EEG reference as zeros (default = [])
%   cfg.refchannel    = cell-array with new EEG reference channel(s), this can be 'all' for a common average reference
%
% See also FT_PREPROCESSING

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
ft_preamble loadvar    data
ft_preamble provenance data

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% the data can be passed as input argument or can be read from disk
hasdata = exist('data', 'var') && ~isempty(data);

% basic check/initialization of input arguments
if ~hasdata
  data = struct([]);
else
  data = ft_checkdata(data);
end

% do a sanity check for incompatible options which are used in ft_preprocessing
cfg = ft_checkconfig(cfg, 'forbidden', {'refmethod', 'montage'});

% set default configuration options
cfg.reref          = ft_getopt(cfg, 'reref', 'yes');
cfg.channel      = ft_getopt(cfg, 'channel', 'all');
cfg.implicitref  = ft_getopt(cfg, 'implicitref');
cfg.refchannel   = ft_getopt(cfg, 'refchannel');

assert(istrue(cfg.reref), 'cannot create a montage without cfg.reref=''yes''');

% here the actual work starts
if hasdata
  cfg.channel = ft_channelselection(cfg.channel, data.label);
  cfg.refchannel = ft_channelselection(cfg.refchannel, cat(1, data.label(:), cfg.implicitref));
else
  if ischar(cfg.channel)
    cfg.channel = {cfg.channel};
  end
  if ischar(cfg.refchannel)
    cfg.refchannel = {cfg.refchannel};
  end
end

% the first montage serves to add the implicit reference channel
montage1 = [];
montage1.labelold = cfg.channel(:);
montage1.labelnew = cfg.channel(:);
if ~isempty(cfg.implicitref)
  montage1.labelnew{end+1} = cfg.implicitref;
end
% the last row for the implicitref will be all zero
montage1.tra = eye(numel(montage1.labelnew), numel(montage1.labelold));

% the second montage serves to subtract the selected reference channels
montage2 = [];
montage2.labelold = montage1.labelnew;
montage2.labelnew = montage1.labelnew;
montage2.tra      = eye(numel(montage2.labelnew));
refsel = match_str(montage2.labelold, cfg.refchannel);
% subtract the mean of the reference channels
montage2.tra(:,refsel) = montage2.tra(:,refsel) - 1/numel(refsel);

% apply montage2 to montage1, the result is the combination of both
montage = ft_apply_montage(montage1, montage2);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble provenance
ft_postamble previous data
ft_postamble history montage

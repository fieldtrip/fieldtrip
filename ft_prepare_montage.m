function [montage, cfg] = ft_prepare_montage(cfg, data)

% FT_PREPARE_MONTAGE creates a referencing scheme based on the input configuration
% options and the channels in the data structure. The resulting montage can be
% given as input to FT_APPLY_MONTAGE, or as cfg.montage to FT_PREPROCESSING.
%
% Use as
%   montage = ft_prepare_montage(cfg, data)
%
% The configuration can contain the following fields:
%   cfg.refmethod     = 'avg', 'bioloar', 'comp' (default = 'avg')
%   cfg.implicitref   = string with the label of the implicit reference, or empty (default = [])
%   cfg.refchannel    = cell-array with new EEG reference channel(s), this can be 'all' for a common average reference
%
% The implicitref option allows adding the implicit reference channel to the data as
% a channel with zeros.
%
% See also FT_PREPROCESSING, FT_APPLY_MONTAGE

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
  data = ft_checkdata(data, 'datatype', {'raw', 'raw+comp', 'timelock', 'timelock+comp'});
end

% do a sanity check for incompatible options which are used in ft_preprocessing and elsewhere
cfg = ft_checkconfig(cfg, 'forbidden', {'montage', 'method'});

% set default configuration options
cfg.refmethod    = ft_getopt(cfg, 'refmethod', 'avg');
cfg.channel      = ft_getopt(cfg, 'channel', 'all');
cfg.implicitref  = ft_getopt(cfg, 'implicitref');
cfg.refchannel   = ft_getopt(cfg, 'refchannel');

if isfield(cfg, 'reref')
  assert(istrue(cfg.reref), 'cannot create a montage without cfg.reref=''yes''');
end

% here the actual work starts
if hasdata
  cfg.channel    = ft_channelselection(cfg.channel, data.label);
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

switch cfg.refmethod
  case 'avg'
    % make a montage that subtracts the mean of all channels specified in cfg.refchannel
    montage2 = [];
    montage2.labelold = montage1.labelnew;
    montage2.labelnew = montage1.labelnew;
    montage2.tra      = eye(numel(montage2.labelnew));
    refsel = match_str(montage2.labelold, cfg.refchannel);
    montage2.tra(:,refsel) = montage2.tra(:,refsel) - 1/numel(refsel);
    
    % apply montage2 to montage1, the result is the combination of both
    montage = ft_apply_montage(montage1, montage2);
    
  case 'bipolar'
    % make a montage for the bipolar derivation of sequential channels
    montage2          = [];
    montage2.labelold = montage1.labelnew;
    montage2.labelnew = strcat(montage1.labelnew(1:end-1),'-',montage1.labelnew(2:end));
    tra_neg           = diag(-ones(numel(montage1.labelnew)-1,1), 1);
    tra_plus          = diag( ones(numel(montage1.labelnew)-1,1),-1);
    montage2.tra      = tra_neg(1:end-1,:)+tra_plus(2:end,:);
    
    % apply montage2 to montage1, the result is the combination of both
    montage = ft_apply_montage(montage1, montage2);
    
  case 'comp'
    assert(isempty(cfg.implicitref), 'implicitref not supported with refmethod=''%s''', cfg.refmethod);
    % construct an linear projection from channels to components
    montage = [];
    montage.labelold = data.topolabel;
    montage.labelnew = data.label;
    montage.tra      = data.unmixing;
    
  case 'invcomp'
    assert(isempty(cfg.implicitref), 'implicitref not supported with refmethod=''%s''', cfg.refmethod);
    % construct an linear projection from components to channels
    montage = [];
    montage.labelold = data.label;
    montage.labelnew = data.topolabel;
    montage.tra      = data.topo;
    
  otherwise
    error('unsupported refmethod=''%s''', cfg.refmethod);
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble provenance
ft_postamble previous data
ft_postamble history montage

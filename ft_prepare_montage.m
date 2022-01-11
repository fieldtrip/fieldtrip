function [montage, cfg] = ft_prepare_montage(cfg, data)

% FT_PREPARE_MONTAGE creates a referencing scheme based on the input configuration
% options and the channels in the data structure. The resulting montage can be
% given as input to FT_APPLY_MONTAGE, or as cfg.montage to FT_PREPROCESSING.
%
% Use as
%   montage = ft_prepare_montage(cfg, data)
%
% The configuration can contain the following fields:
%   cfg.refmethod     = 'avg', 'comp', 'bipolar', 'laplace', 'doublebanana', 'longitudinal', 'circumferential', 'transverse' (default = 'avg')
%   cfg.implicitref   = string with the label of the implicit reference, or empty (default = [])
%   cfg.refchannel    = cell-array with new EEG reference channel(s), this can be 'all' for a common average reference
%   cfg.groupchans    = 'yes' or 'no', should channels be rereferenced in separate groups
%                       for bipolar and laplace methods, this requires channnels to be
%                       named using an alphanumeric code, where letters represent the
%                       group and numbers represent the order of the channel whithin
%                       its group (default = 'no')
%
% The implicitref option allows adding the implicit reference channel to the data as
% a channel with zeros.
%
% The resulting montage is a structure with the fields
%   montage.tra      = MxN matrix
%   montage.labelold = Nx1 cell-array
%   montage.labelnew = Mx1 cell-array
%
% As an example, an output bipolar montage could look like this
%   bipolar.labelold  = {'1',   '2',   '3',   '4'}
%   bipolar.labelnew  = {'1-2', '2-3', '3-4'}
%   bipolar.tra       = [
%     +1 -1  0  0
%      0 +1 -1  0
%      0  0 +1 -1
%   ];
%
% See also FT_PREPROCESSING, FT_APPLY_MONTAGE

% Copyright (C) 2017-2022, Robert Oostenveld
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
cfg.groupchans   = ft_getopt(cfg, 'groupchans', 'no');

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
% the last row for the implicitref will be all zero
montage1 = [];
montage1.labelold = cfg.channel(:);
montage1.labelnew = cfg.channel(:);
if ~isempty(cfg.implicitref)
  montage1.labelnew{end+1} = cfg.implicitref;
end
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
    
  case 'bipolar'
    if istrue(cfg.groupchans)
      % this uses recursion to make the montage for each group
      label = montage1.labelnew;
      group = cell(size(label));
      for i=1:numel(label)
        group(i) = regexp(label{i}, '^[a-zA-Z]+', 'match');
      end
      group = unique(group, 'stable');
      
      tmpcfg = keepfields(cfg, {'refmethod', 'implicitref', 'refchannel', 'channel'});
      tmpcfg.groupchans = 'no';
      tmpcfg.showcallinfo = 'no';
      
      groupmontage = cell(size(group));
      for g = 1:length(group) % for each group of channels
        % select the channels within the group
        tmpcfg.channel = label(startsWith(label, group{g}));
        groupmontage{g} = ft_prepare_montage(tmpcfg);
      end
      
      % start with a montage that keeps all channels the same
      montage2 = [];
      montage2.labelold = label(:)';
      montage2.labelnew = label(:);
      montage2.tra = eye(length(label));
      
      for g = length(group):-1:1
        % apply the montage for each groep (i.e., combine them), keep all other channels as they are
        % unused channels end up at the end, hence doing it in reverse order keeps the channel order consistent
        montage2 = ft_apply_montage(montage2, groupmontage{g}, 'keepunused', true);
      end
      
    else
      % create a montage for the bipolar derivation of sequential channels
      montage2          = [];
      montage2.labelold = montage1.labelnew;
      montage2.labelnew = strcat(montage1.labelnew(1:end-1),'-',montage1.labelnew(2:end));
      tra_neg           = diag(-ones(numel(montage1.labelnew)-1,1), 1);
      tra_plus          = diag( ones(numel(montage1.labelnew)-1,1),-1);
      montage2.tra      = tra_neg(1:end-1,:)+tra_plus(2:end,:);
      
    end % if groupchans
    
    % apply montage2 to montage1, the result is the combination of both
    montage = ft_apply_montage(montage1, montage2);
    
  case 'laplace'
    if istrue(cfg.groupchans)
      % this uses recursion to make the montage for each group
      label = montage1.labelnew;
      group = cell(size(label));
      for i=1:numel(label)
        group(i) = regexp(label{i}, '^[a-zA-Z]+', 'match');
      end
      group = unique(group, 'stable');
      
      tmpcfg = keepfields(cfg, {'refmethod', 'implicitref', 'refchannel', 'channel'});
      tmpcfg.groupchans = 'no';
      tmpcfg.showcallinfo = 'no';
      
      groupmontage = cell(size(group));
      for g = 1:length(group) % for each group of channels
        % select the channels within the group
        tmpcfg.channel = label(startsWith(label, group{g}));
        groupmontage{g} = ft_prepare_montage(tmpcfg);
      end
      
      % start with a montage that keeps all channels the same
      montage2 = [];
      montage2.labelold = label(:)';
      montage2.labelnew = label(:);
      montage2.tra = eye(length(label));
      
      for g = length(group):-1:1
        % apply the montage for each groep (i.e., combine them), keep all other channels as they are
        % unused channels end up at the end, hence doing it in reverse order keeps the channel order consistent
        montage2 = ft_apply_montage(montage2, groupmontage{g}, 'keepunused', true);
      end
      
    else
      % create a montage for the laplacian derivation of neighboring channels where
      % each channel is re-referenced against the mean of the two nearest channels,
      % while channels at the extremities will be re-referenced against their closest
      % neighbor; this is particularly useful for sEEG shafts
      
      montage2          = [];
      montage2.labelold = montage1.labelnew;
      tra_neg1          = diag(-0.5*ones(numel(montage1.labelnew)-1,1), 1);
      tra_neg2          = diag(-0.5*ones(numel(montage1.labelnew)-1,1), -1);
      
      if numel(montage2.labelold) > 2
        montage2.labelnew = montage1.labelnew;
        tra_plus          = diag(ones(numel(montage1.labelnew),1),0);
        montage2.tra      = tra_neg1+tra_neg2+tra_plus;
        montage2.tra(1,2) = -1;
        montage2.tra(end,end-1) = -1;
      else
        % The following commands are just meant to avoid errors, in case
        % only one channel is given as input: the algorithm will behave
        % exactly as it would in the 'bipolar' case
        montage2.labelnew = cell(numel(montage1.labelnew),0);
        tra_plus          = diag(ones(numel(montage1.labelnew)-1,1),-1);
        montage2.tra      = tra_neg1+tra_neg2+tra_plus(2:end,:);
      end
      
    end % if groupchans
    
    % apply montage2 to montage1, the result is the combination of both
    montage = ft_apply_montage(montage1, montage2);
    
  case {'longitudinal', 'doublebanana'}
    % see https://www.learningeeg.com/montages-and-technical-components
    montage = [];
    montage.labelold = {'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'T3', 'C3', 'Cz', 'C4', 'T4', 'T5', 'P3', 'Pz', 'P4', 'T6', 'O1', 'Oz', 'O2'};
    montage.labelnew = {
      % left temporal chain
      'Fp1-F7'
      'F7-T3'
      'T3-T5'
      'T5-O1'
      % left parasagittal chain
      'Fp1-F3'
      'F3-C3'
      'C3-P3'
      'P3-O1'
      % central chain
      'Fz-Cz'
      'Cz-Pz'
      % right parasagittal chain
      'Fp2-F4'
      'F4-C4'
      'C4-P4'
      'P4-O2'
      % right temporal chain
      'Fp2-F8'
      'F8-T4'
      'T4-T6'
      'T6-O2'
      };
    % construct the montage from new channel labels
    montage.tra = zeros(length(montage.labelnew), length(montage.labelold));
    for i=1:length(montage.labelnew)
      lab = split(montage.labelnew{i}, '-');
      montage.tra(i, strcmp(montage.labelold, lab{1})) = +1;
      montage.tra(i, strcmp(montage.labelold, lab{2})) = -1;
    end
    % do a sanity check on the montage
    assert(all(sum(montage.tra, 2)==0));
    assert(all(max(montage.tra, [], 2)==+1));
    assert(all(min(montage.tra, [], 2)==-1));
    
  case 'circumferential'
    % see https://www.learningeeg.com/montages-and-technical-components
    montage = [];
    montage.labelold = {'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'T3', 'C3', 'Cz', 'C4', 'T4', 'T5', 'P3', 'Pz', 'P4', 'T6', 'O1', 'Oz', 'O2'};
    montage.labelnew = {
      'Fp1-F7'
      'F7-T3'
      'T3-T5'
      'T5-O1'
      'O1-O2'
      'O2-T6'
      'T6-T4'
      'T4-F8'
      'F8-Fp2'
      'Fp2-Fp1'
      };
    % construct the montage from new channel labels
    montage.tra = zeros(length(montage.labelnew), length(montage.labelold));
    for i=1:length(montage.labelnew)
      lab = split(montage.labelnew{i}, '-');
      montage.tra(i, strcmp(montage.labelold, lab{1})) = +1;
      montage.tra(i, strcmp(montage.labelold, lab{2})) = -1;
    end
    % do a sanity check on the montage
    assert(all(sum(montage.tra, 2)==0));
    assert(all(max(montage.tra, [], 2)==+1));
    assert(all(min(montage.tra, [], 2)==-1));
    
  case 'transverse'
    % this is inspired by https://doi.org/10.1016/j.earlhumdev.2011.08.008 but probably better documented elsewhere
    % see https://www.learningeeg.com/montages-and-technical-components
    % this particular implementation only considers the 10% distances
    montage = [];
    montage.labelold = {'Fp1', 'Fpz', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'T3', 'C3', 'Cz', 'C4', 'T4', 'T5', 'P3', 'Pz', 'P4', 'T6', 'O1', 'Oz', 'O2'};
    montage.labelnew = {
      'Fp1-Fp2'
      'F7-F3'
      'F3-Fz'
      'Fz-F4'
      'F4-F8'
      'T3-C3'
      'C3-Cz'
      'Cz-C4'
      'C4-T4'
      'T5-P3'
      'P3-Pz'
      'Pz-P4'
      'P4-T6'
      'O1-O2'
      };
    % construct the montage from new channel labels
    montage.tra = zeros(length(montage.labelnew), length(montage.labelold));
    for i=1:length(montage.labelnew)
      lab = split(montage.labelnew{i}, '-');
      montage.tra(i, strcmp(montage.labelold, lab{1})) = +1;
      montage.tra(i, strcmp(montage.labelold, lab{2})) = -1;
    end
    % do a sanity check on the montage
    assert(all(sum(montage.tra, 2)==0));
    assert(all(max(montage.tra, [], 2)==+1));
    assert(all(min(montage.tra, [], 2)==-1));
    
  otherwise
    error('unsupported refmethod=''%s''', cfg.refmethod);
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble provenance
ft_postamble previous data
ft_postamble history montage

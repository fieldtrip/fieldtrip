function [obj] = convert_event(obj, target, varargin)

% CONVERT_EVENT converts between different representations of events, trials or
% artifacts.
%
% Use as
%   [object] = convert_event(object, target, ....)
%
% The different representations that can be used for the input object are
% 1) event structure, see FT_READ_EVENT
% 2) matrix representation as in trl (Nx3), see FT_DEFINETRIAL
% 3) matrix representation as in artifact (Nx2), see FT_ARTIFACT_xxx
% 4) continuously sampled boolean vector representation that represents the onset and offset of the events or artifact
%
% The different representations that can be specified for the target representation
% are 'event', 'trl', 'artifact', 'boolvec'
%
% Additional options should be specified in key-value pairs and can be
%   'endsample'  = scalar
%   'typenames'  = cell-array with strings
%   'valuenames' = cell-array with strings
%
% See FT_READ_EVENT, FT_DEFINETRIAL, FT_REJECTARTIFACT, FT_ARTIFACT_xxx

% Copyright (C) 2009, Ingrid Nieuwenhuis
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

% Check if target is specified correctly
if ~ismember({'event', 'trl', 'artifact', 'boolvec'}, target)
  ft_error('target has to be ''event'', ''trl'', ''artifact'', or ''boolvec''.')
end

% Get the options
endsample  = ft_getopt(varargin, 'endsample');
typenames  = ft_getopt(varargin, 'typenames');
valuenames = ft_getopt(varargin, 'valuenames');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine what the input object is
if isempty(obj)
  input = 'empty';
  
elseif isstruct(obj)
  input = 'event';
  
elseif iscell(obj)
  if isempty(obj{1})
    input = 'empty';
  elseif size(obj{1},2) == 3
    input = 'trl';
  elseif size(obj{1},2) == 2
    input = 'artifact';
  elseif size(obj{1},2) > 3
    % could be a strange trl-matrix with multiple columns
    input = 'trl';
    for i = 1:length(obj)
      if ~isempty(obj{i})
        obj{i} = obj{i}(:,1:3);
      end
    end
  else
    ft_error('incorrect input object, see help for what is allowed.')
  end
  
elseif islogical(obj)
  input = 'boolvec';
  
elseif size(obj,2) == 3
  input = 'trl';
  
elseif size(obj,2) == 2
  input = 'artifact';
  
elseif size(obj,2) > 3
  tmp = unique(obj);
  if isempty(find(tmp>2, 1))
    input = 'boolvec';
    obj = logical(obj);
  else
    % it is at least not boolean, but it could be a strange trl-matrix with multiple columns
    input = 'trl';
    obj = obj(:,1:3);
  end
  
else
  ft_error('unrecognized input object')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the conversion
if (strcmp(input, 'trl') || strcmp(input, 'artifact') || strcmp(input, 'empty')) && strcmp(target, 'boolvec')
  if ~isempty(endsample)
    obj = artifact2boolvec(obj,endsample);
  else
    obj = artifact2boolvec(obj);
  end
  
elseif strcmp(input, 'boolvec') && strcmp(target,'artifact' )
  obj = boolvec2artifact(obj);
  
elseif strcmp(input, 'boolvec') && strcmp(target,'trl' )
  obj = boolvec2artifact(obj);
  if iscell(obj)
    for i=1:length(obj)
      obj{i}(:,3) = 0;
    end
  else
    obj(:,3) = 0;
  end
  
elseif (strcmp(input, 'trl') || strcmp(input, 'artifact') || strcmp(input, 'empty')) && strcmp(target, 'event')
  obj = artifact2event(obj, typenames);
  
elseif strcmp(input, 'event') && strcmp(target, 'boolvec')
  obj = event2boolvec(obj, typenames, valuenames, endsample);
  
elseif strcmp(input, 'artifact') && strcmp(target, 'trl')
  % add an extra column with the offset
  if iscell(obj)
    for i=1:length(obj)
      obj{i}(:,3) = 0;
    end
  else
    obj(:,3) = 0;
  end
  
elseif strcmp(input, 'trl') && strcmp(target, 'artifact')
  % remove the extra column with the offset
  if iscell(obj)
    for i=1:length(obj)
      obj{i}(:,3:end) = [];
    end
  else
    obj(:,3:end) = [];
  end
  
elseif strcmp(input, 'empty')
  obj = [];
  
else
  ft_warning('conversion from %s to %s is not supported yet', input, target);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function makes a Boolean vector out of an event structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function boolvec = event2boolvec(event, typenames, valuenames, endsample)

if nargin<3 || isempty(endsample)
  if isempty([event.duration])
    endsample = [event.sample];
  elseif all([event.duration]==0)
    % the duration is infinitely short
    endsample = [event.sample];
  elseif all([event.duration]==1)
    % the duration is one sample
    endsample = [event.sample];
  else
    endsample = [event.sample] + [event.duration] - 1;
  end
  endsample = max(endsample);
end

boolvec = false(length(typenames), endsample);

for i=1:length(typenames)
  sel = strcmp({event.type}, typenames{i});
  begsample = [event(sel).sample];
  duration  = [event(sel).duration];
  if isempty(duration)
    endsample = begsample;
  elseif all(duration==0)
    endsample = begsample;
  elseif all(duration==1)
    endsample = begsample;
  else
    endsample = begsample + duration - 1;
  end
  for j=1:length(begsample)
    boolvec(i, begsample(j):endsample(j)) = true;
  end
end

for i=1:length(valuenames)
  sel = strcmp({event.value}, valuenames{i});
  begsample = [event(sel).sample];
  duration  = [event(sel).duration];
  if isempty(duration)
    endsample = begsample;
  elseif all(duration==0)
    endsample = begsample;
  elseif all(duration==1)
    endsample = begsample;
  else
    endsample = begsample + duration - 1;
  end
  for j=1:length(begsample)
    boolvec(i, begsample(j):endsample(j)) = true;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function makes a Boolean vector (or matrix when artifact is cell-array of
% multiple artifact definitions) with 0 for artifact free sample and 1 for sample
% containing an artifact according to artifact specification. Length of vector is
% (from sample 1) last sample as defined in the artifact definition, or when
% datendsample is specified vector is length datendsample.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function boolvec = artifact2boolvec(artifact, endsample)

if nargin==1
  if ~iscell(artifact) % assume only one artifact is given
    if isempty(artifact)
      ft_error('When input object is empty ''endsample'' must be specified to convert into boolvec')
    else
      endsample = max(artifact(:,2));
    end
  elseif length(artifact) == 1
    if isempty(artifact{1})
      ft_error('When input object is empty ''endsample'' must be specified to convert into boolvec')
    else
      endsample = max(artifact{1}(:,2));
    end
  else
    ft_error('when giving multiple artifact definitions, endsample should be specified to assure all output vectors are of the same length')
  end
end

if ~iscell(artifact)
  artifact = {artifact};
end

% make boolvec
boolvec = zeros(length(artifact), endsample);
breakflag = 0;
for i=1:length(artifact)
  for j=1:size(artifact{i},1)
    artbegsample = artifact{i}(j,1);
    artendsample = artifact{i}(j,2);
    if artbegsample > endsample
      ft_warning('artifact definition contains later samples than endsample, these samples are ignored')
      break
    elseif artendsample > endsample
      ft_warning('artifact definition contains later samples than endsample, these samples are ignored')
      artendsample = endsample;
      breakflag = 1;
    end
    boolvec(i, artbegsample:artendsample) = 1;
    if breakflag
      break
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function makes an artifact definition or a cell-array of artifact
% definitions from a Boolean vector (or matrix).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function artifact = boolvec2artifact(boolvec)

for i=1:size(boolvec,1)
  tmp = diff([0 boolvec(i,:) 0]);
  artbeg = find(tmp==+1);
  artend = find(tmp==-1) - 1;
  artifact{i} = [artbeg' artend'];
end

if length(artifact) == 1
  artifact = artifact{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function makes an event structure from artifact definition (or cell array of
% artifact definitions). event.type is always 'artifact', but in case of cell-arrays
% of artifacts is is also possible to hand 'typenames' with length(artifact)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function event = artifact2event(artifact, typenames)

if ~iscell(artifact)
  artifact = {artifact};
end

if ~isempty(typenames)
  if length(artifact) ~= length(typenames)
    ft_error('length typenames should be the same as length artifact')
  end
end

event = [];
for i=1:length(artifact)
  for j=1:size(artifact{i},1)
    event(end+1).sample   = artifact{i}(j,1);
    event(end  ).duration = artifact{i}(j,2)-artifact{i}(j,1)+1;
    if ~isempty(typenames)
      event(end).type     = typenames{i};
    elseif size(artifact{i},2) == 2
      event(end).type     = 'artifact';
    elseif size(artifact{i},2) == 3
      event(end).type     = 'trial';
    end
    event(end  ).value    = [];
    event(end  ).offset   = [];
  end
end

function [obj] = convert_event(obj, target, varargin)

% CONVERT_EVENT converts between the different representations of events, 
% which can be
%   1) event structure, see FT_READ_EVENT
%   2) matrix representation as in trl (Nx3), see FT_DEFINETRIAL
%   3) matrix representation as in artifact (Nx2), see FT_ARTIFACT_xxx
%   4) boolean vector representation with 1 for samples containing trial/artifact
%
% Use as
%   [object] = convert_event(object, target, ....)
%
% Possible input objects types are
%   event structure
%   Nx3 trl matrix, or cell array of multiple trl definitions
%   Nx2 artifact matrix, or cell array of multiple artifact definitions
%   boolean vector, or matrix with multiple vectors as rows
%
% Possible targets are 'event', 'trl', 'artifact', 'boolvec'
%
% Additional options should be specified in key-value pairs and can be
%   'endsample' = 
%   'typenames' = 
%
% See READ_EVENT, DEFINETRIAL, REJECTARTIFACT, ARTIFACT_xxx

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
if sum(strcmp(target, {'event', 'trl', 'artifact', 'boolvec'})) < 1
  error('target has to be ''event'', ''trl'', ''artifact'', or ''boolvec''.')
end

% Get the options
endsample = ft_getopt(varargin, 'endsample');
typenames = ft_getopt(varargin, 'typenames');

% Determine what the input object is
if isempty(obj)
  input_obj = 'empty';
elseif isstruct(obj)
  input_obj = 'event';
elseif iscell(obj)
  if isempty(obj{1})
    input_obj = 'empty';
  elseif size(obj{1},2) == 3
    input_obj = 'trl';
  elseif size(obj{1},2) == 2
    input_obj = 'artifact';
  elseif size(obj{1},2) > 3
    % could be a strange trl-matrix with multiple columns
    input_obj = 'trl';
    for i = 1:length(obj)
      if ~isempty(obj{i})
        obj{i} = obj{i}(:,1:3);
      end
    end
  else
    error('incorrect input object, see help for what is allowed.')
  end
elseif islogical(obj)
  input_obj = 'boolvec';
elseif size(obj,2) == 3
  input_obj = 'trl';
elseif size(obj,2) == 2
  input_obj = 'artifact';
elseif size(obj,2) > 3
  tmp = unique(obj);
  if isempty(find(tmp>2, 1))
    input_obj = 'boolvec';
    obj = logical(obj);
  else
    %it is at least not boolean but could be a strange
    %trl-matrix with multiple columns
    input_obj = 'trl';
    obj = obj(:,1:3);
  end
else
  error('incorrect input object, see help for what is allowed.')
end

% do conversion
if (strcmp(input_obj, 'trl') || strcmp(input_obj, 'artifact') || strcmp(input_obj, 'empty')) && strcmp(target, 'boolvec')
  if ~isempty(endsample)
    obj = artifact2artvec(obj,endsample);
  else
    obj = artifact2artvec(obj);
  end
elseif strcmp(input_obj, 'boolvec') && strcmp(target,'artifact' )
  obj = artvec2artifact(obj);
elseif strcmp(input_obj, 'boolvec') && strcmp(target,'trl' )
  obj = artvec2artifact(obj);
  if iscell(obj)
    for i=1:length(obj)
      obj{i}(:,3) = 0;
    end
  else
    obj(:,3) = 0;
  end
elseif (strcmp(input_obj, 'trl') || strcmp(input_obj, 'artifact')) && strcmp(target, 'event')
  obj = artifact2event(obj, typenames);
elseif strcmp(input_obj, 'artifact') && strcmp(target,'trl')
  if iscell(obj)
    for i=1:length(obj)
      obj{i}(:,3) = 0;
    end
  else
    obj(:,3) = 0;
  end
elseif strcmp(input_obj, 'trl') && strcmp(target,'artifact')
  if iscell(obj)
    for i=1:length(obj)
      obj{i}(:,3:end) = [];
    end
  else
    obj(:,3:end) = [];
  end
elseif strcmp(input_obj, 'empty')
  obj = [];
else
  warning('conversion not supported yet') %FIXME
end


%%%%%%%%%%%%%%% SUBFUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function artvec = artifact2artvec(varargin)
% ARTIFACT2ARTVEC makes boolian vector (or matrix when artifact is
% cell array of multiple artifact definitions) with 0 for artifact free
% sample and 1 for sample containing an artifact according to artifact
% specification. Length of vector is (from sample 1) last sample as
% defined in the artifact definition, or when datendsample is speciefied
% vector is length datendsample.

artifact = varargin{1};
if length(varargin) == 1
  if ~iscell(artifact) % assume only one artifact is given
    if isempty(artifact)
      error('When input object is empty ''endsample'' must be specified to convert into boolvec')
    else
      endsample = max(artifact(:,2));
    end
  elseif length(artifact) == 1
    if isempty(artifact{1})
      error('When input object is empty ''endsample'' must be specified to convert into boolvec')
    else
      endsample = max(artifact{1}(:,2));
    end
  else
    error('when giving multiple artifact definitions, endsample should be specified to assure all output vectors are of the same length')
  end
elseif length(varargin) == 2
  endsample = varargin{2};
elseif length(varargin) > 2
  error('too many input arguments')
end
if ~iscell(artifact)
  artifact = {artifact};
end

% make artvec
artvec = zeros(length(artifact), endsample);
breakflag = 0;
for i=1:length(artifact)
  for j=1:size(artifact{i},1)
    artbegsample = artifact{i}(j,1);
    artendsample = artifact{i}(j,2);
    if artbegsample > endsample
      warning('artifact definition contains later samples than endsample, these samples are ignored')
      break
    elseif artendsample > endsample
      warning('artifact definition contains later samples than endsample, these samples are ignored')
      artendsample = endsample;
      breakflag = 1;
    end
    artvec(i, artbegsample:artendsample) = 1;
    if breakflag
      break
    end
  end
end

%%%%%%%%%%%%%%% SUBFUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function artifact = artvec2artifact(artvec)
% ARTVEC2ARTIFACT makes artifact definition (or cell array of artifact
% definitions) from boolian vector (or matrix) with [artbegsample
% artendsample]. Assumed is that the artvec starts with sample 1.

for i=1:size(artvec,1)
  tmp = diff([0 artvec(i,:) 0]);
  artbeg = find(tmp==+1);
  artend = find(tmp==-1) - 1;
  artifact{i} = [artbeg' artend'];
end

if length(artifact) == 1
  artifact = artifact{1};
end

%%%%%%%%%%%%%%% SUBFUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function event = artifact2event(artifact, typenames)
% ARTIFACT2EVENT makes event structure from artifact definition (or cell
% array of artifact definitions). event.type is always 'artifact', but
% incase of cellarays of artifacts is is also possible to hand 'typenames'
% with length(artifact)

if ~iscell(artifact)
  artifact = {artifact};
end

if ~isempty(typenames)
  if length(artifact) ~= length(typenames)
    error('length typenames should be the same as length artifact')
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



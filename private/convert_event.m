function [obj] = convert_event(obj, target, varargin)

% CONVERT_EVENT converts between the different representations of events.
% The representations are:
%   event structure (see READ_EVENT)
%   matrix representation as in trl (Nx3) or artifact (Nx2), with
%   [begsample endsample] (see DEFINETRIAL, ARTIFACT_xxx functions)
%   boolean vector representation with 1 for samples containing trial/artifact
%
% Use as
%   [object] = convert_event(object, target, ....)
%
% The following input objects are supported:
%   event structure
%   trl (or cell array of multiple trl definitions)
%   artifact (or cell array of multiple artifact definitions)
%   boolian vector (or matrix of multiple vectors)
%
% Possible targets are 'event', 'trl', 'artifact', 'boolvec'
%
% Additional options should be specified in key-value pairs and can be
%   'endsample'
%
% See READ_EVENT, DEFINETRIAL, REJECTARTIFACT, ARTIFACT_xxx

% Copyright (C) 2009, Ingrid Nieuwenhuis
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% Check if target is specified correctly
if sum(strcmp(target, {'event', 'trl', 'artifact', 'boolvec'})) < 1
  error('target has to be ''event'', ''trl'', ''artifact'', or ''boolvec''.')
end

% Get the options
endsample = keyval('endsample',  varargin);

% Determine what the input object is
if isstruct(obj)
  input_obj = 'event';
elseif iscell(obj)
  if size(obj{1},2) == 3
    input_obj = 'trl';
  elseif size(obj{1},2) == 2
    input_obj = 'artifact';
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
  if length(tmp) == 2 && tmp(1) == 0 && tmp(2) == 1
    input_obj = 'boolvec';
    obj = logical(obj);
  end
else
  error('incorrect input object, see help for what is allowed.')
end

% do conversion
if (strcmp(input_obj, 'trl') || strcmp(input_obj, 'artifact')) && strcmp(target, 'boolvec')
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
    endsample = max(artifact(:,2));
  elseif length(artifact) == 1
    endsample = max(artifact{1}(:,2));
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


function [sens] = apply_montage(sens, montage, varargin)

% APPLY_MONTAGE changes the montage of an electrode or gradiometer array. A
% montage can be used for EEG rereferencing, MEG synthetic gradients, MEG
% planar gradients or unmixing using ICA. This function applies the montage
% to the sensor array. The sensor array can subsequently be used for
% forward computation and source reconstruction of the data.
%
% Use as
%   [sens]    = apply_montage(sens,     montage,  ...)
%   [data]    = apply_montage(data,     montage,  ...)
%   [montage] = apply_montage(montage1, montage2, ...)
% where the input is a FieldTrip sensor definition as obtained from READ_SENS
% or a FieldTrip raw data structure as obtained from PREPROCESSING.
%
% A montage is specified as a structure with the fields
%   montage.tra      = MxN matrix
%   montage.labelnew = Mx1 cell-array
%   montage.labelorg = Nx1 cell-array
%
% Additional options should be specified in key-value pairs and can be
%   'keepunused'    string, 'yes' or 'no' (default = 'no')
%   'inverse'       string, 'yes' or 'no' (default = 'no')
%
% If the first input is a montage, then the second input montage will be
% applied to the first. In effect the resulting montage will first do
% montage1, then montage2.
%
% See also READ_SENS, TRANSFORM_SENS

% Copyright (C) 2008, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% get optional input arguments
keepunused    = keyval('keepunused',    varargin{:}); if isempty(keepunused),    keepunused    = 'no';  end
inverse       = keyval('inverse',       varargin{:}); if isempty(inverse),       inverse       = 'no';  end

% check the consistency of the input sensor array or data
if isfield(sens, 'labelorg') && isfield(sens, 'labelnew')
  % the input data structure is also a montage, i.e. apply the montages sequentially
  sens.label = sens.labelnew;
end

% check the consistency of the montage
if size(montage.tra,1)~=length(montage.labelnew)
  error('the number of channels in the montage is inconsistent');
elseif size(montage.tra,2)~=length(montage.labelorg)
  error('the number of channels in the montage is inconsistent');
end

if strcmp(inverse, 'yes')
  % apply the inverse montage, i.e. undo a previously applied montage
  tmp.labelnew = montage.labelorg; % swap around
  tmp.labelorg = montage.labelnew; % swap around
  tmp.tra      = full(montage.tra);
  if rank(tmp.tra) < length(tmp.tra)
    warning('the linear projection for the montage is not full-rank, the resulting data will have reduced dimensionality');
    tmp.tra      = pinv(tmp.tra);
  else
    tmp.tra      = inv(tmp.tra);
  end
  montage      = tmp;
end

% use default transfer from sensors to channels if not specified
if ~isfield(sens, 'tra') && ~isfield(sens, 'trial')
  nchan = size(sens.pnt,1);
  sens.tra = sparse(eye(nchan));
end

% select and keep the columns that are non-empty, i.e. remove the empty columns
selcol           = find(~all(montage.tra==0, 1));
montage.tra      = montage.tra(:,selcol);
montage.labelorg = montage.labelorg(selcol);
clear selcol

% select and remove the columns corresponding to channels that are not present in the original data
remove = setdiff(montage.labelorg, intersect(montage.labelorg, sens.label));
selcol = match_str(montage.labelorg, remove);
% we cannot just remove the colums, all rows that depend on it should also be removed
selrow = false(length(montage.labelnew),1);
for i=1:length(selcol)
  selrow = selrow & (montage.tra(:,selcol(i))~=0);
end
% convert from indices to logical vector
selcol = indx2logical(selcol, length(montage.labelorg));
% remove rows and columns
montage.labelorg = montage.labelorg(~selcol);
montage.labelnew = montage.labelnew(~selrow);
montage.tra = montage.tra(~selrow, ~selcol);
clear remove selcol selrow i
% add columns for the channels that are present in the data but not involved in the montage, and stick to the original order in the data
[add, ix] = setdiff(sens.label, montage.labelorg);
add = sens.label(sort(ix));
m = size(montage.tra,1);
n = size(montage.tra,2);
k = length(add);
if strcmp(keepunused, 'yes')
  % add the channels that are not rereferenced to the input and output
  montage.tra((m+(1:k)),(n+(1:k))) = eye(k);
  montage.labelorg = cat(1, montage.labelorg(:), add(:));
  montage.labelnew = cat(1, montage.labelnew(:), add(:));
else
  % add the channels that are not rereferenced to the input montage only
  montage.tra(:,(n+(1:k))) = zeros(m,k);
  montage.labelorg = cat(1, montage.labelorg(:), add(:));
end
clear add m n k

% determine whether all channels are unique
m = size(montage.tra,1);
n = size(montage.tra,2);
if length(unique(montage.labelnew))~=m
  error('not all output channels of the montage are unique');
end
if length(unique(montage.labelorg))~=n
  error('not all input channels of the montage are unique');
end

% determine whether all channels that have to be rereferenced are available
if length(intersect(sens.label, montage.labelorg))~=length(montage.labelorg)
  error('not all channels that are required in the montage are available in the data');
end

% reorder the columns of the montage matrix
[selsens, selmont] = match_str(sens.label, montage.labelorg);
montage.tra        = sparse(montage.tra(:,selmont));
montage.labelorg   = montage.labelorg(selmont);

if isfield(sens, 'labelorg') && isfield(sens, 'labelnew')
  % apply the montage on top of the other montage
  sens       = rmfield(sens, 'label');
  sens.tra   = montage.tra * sens.tra;
  sens.labelnew = montage.labelnew;

elseif isfield(sens, 'tra')
  % apply the montage to the sensor array
  sens.tra   = montage.tra * sens.tra;
  sens.label = montage.labelnew;

elseif isfield(sens, 'trial')
  % apply the montage to the data that was preprocessed using fieldtrip
  Ntrials = numel(sens.trial);
  for i=1:Ntrials
    fprintf('processing trial %d from %d\n', i, Ntrials);
    if isa(sens.trial{i}, 'single')
      % sparse matrices and single precision do not match
      sens.trial{i}   = full(montage.tra) * sens.trial{i};
    else
      sens.trial{i}   = montage.tra * sens.trial{i};
    end
  end
  sens.label = montage.labelnew;
else
  error('unrecognized input');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = indx2logical(x, n)
y = false(1,n);
y(x) = true;

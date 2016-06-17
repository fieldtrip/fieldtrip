function [input] = ft_apply_montage(input, montage, varargin)

% FT_APPLY_MONTAGE changes the montage of an electrode or gradiometer array. A
% montage can be used for EEG rereferencing, MEG synthetic gradients, MEG
% planar gradients or unmixing using ICA. This function applies the montage
% to the input EEG or MEG sensor array, which can subsequently be used for
% forward computation and source reconstruction of the data.
%
% Use as
%   [sens]    = ft_apply_montage(sens,     montage,  ...)
%   [data]    = ft_apply_montage(data,     montage,  ...)
%   [freq]    = ft_apply_montage(freq,     montage,  ...)
%   [montage] = ft_apply_montage(montage1, montage2, ...)
%
% A montage is specified as a structure with the fields
%   montage.tra      = MxN matrix
%   montage.labelold = Nx1 cell-array
%   montage.labelnew = Mx1 cell-array
%
% As an example, a bipolar montage could look like this
%   bipolar.labelold  = {'1', '2', '3', '4'}
%   bipolar.labelnew  = {'1-2', '2-3', '3-4'}
%   bipolar.tra       = [
%     +1 -1  0  0
%      0 +1 -1  0
%      0  0 +1 -1
%   ];
%
% The montage can optionally also specify the channel type and unit of the input
% and output data with
%   montage.chantypeold = Nx1 cell-array
%   montage.chantypenew = Mx1 cell-array
%   montage.chanunitold = Nx1 cell-array
%   montage.chanunitnew = Mx1 cell-array
%
% Additional options should be specified in key-value pairs and can be
%   'keepunused'    string, 'yes' or 'no' (default = 'no')
%   'inverse'       string, 'yes' or 'no' (default = 'no')
%   'balancename'   string, name of the montage (default = '')
%   'feedback'      string, see FT_PROGRESS (default = 'text')
%   'warning'       boolean, whether to show warnings (default = true)
%
% If the first input is a montage, then the second input montage will be
% applied to the first. In effect, the output montage will first do
% montage1, then montage2.
%
% See also FT_READ_SENS, FT_TRANSFORM_SENS

% Copyright (C) 2008-2016, Robert Oostenveld
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

if iscell(input) && iscell(input)
  % this represents combined EEG, ECoG and/or MEG
  for i=1:numel(input)
    input{i} = ft_apply_montage(input{i}, montage, varargin{:});
  end
  return
end

% use "old/new" instead of "org/new"
montage = fixmontage(montage);
input   = fixmontage(input); % the input might also be a montage

% get optional input arguments
keepunused  = ft_getopt(varargin, 'keepunused',  'no');
inverse     = ft_getopt(varargin, 'inverse',     'no');
feedback    = ft_getopt(varargin, 'feedback',    'text');
showwarning = ft_getopt(varargin, 'warning',     true);
bname       = ft_getopt(varargin, 'balancename', '');

if istrue(showwarning)
  warningfun = @warning;
else
  warningfun = @nowarning;
end

% these are optional, at the end we will clean up the output in case they did not exist
haschantype = (isfield(input, 'chantype') || isfield(input, 'chantypenew')) && all(isfield(montage, {'chantypeold', 'chantypenew'}));
haschanunit = (isfield(input, 'chanunit') || isfield(input, 'chanunitnew')) && all(isfield(montage, {'chanunitold', 'chanunitnew'}));

% make sure they always exist to facilitate the remainder of the code
if ~isfield(montage, 'chantypeold')
  montage.chantypeold = repmat({'unknown'}, size(montage.labelold));
  if isfield(input, 'chantype') && ~istrue(inverse)
    warning('copying input chantype to montage');
    [sel1, sel2] = match_str(montage.labelold, input.label);
    montage.chantypeold(sel1) = input.chantype(sel2);
  end
end

if ~isfield(montage, 'chantypenew')
  montage.chantypenew = repmat({'unknown'}, size(montage.labelnew));
  if isfield(input, 'chantype') && istrue(inverse)
    warning('copying input chantype to montage');
    [sel1, sel2] = match_str(montage.labelnew, input.label);
    montage.chantypenew(sel1) = input.chantype(sel2);
  end
end

if ~isfield(montage, 'chanunitold')
  montage.chanunitold = repmat({'unknown'}, size(montage.labelold));
  if isfield(input, 'chanunit') && ~istrue(inverse)
    warning('copying input chanunit to montage');
    [sel1, sel2] = match_str(montage.labelold, input.label);
    montage.chanunitold(sel1) = input.chanunit(sel2);
  end
end

if ~isfield(montage, 'chanunitnew')
  montage.chanunitnew = repmat({'unknown'}, size(montage.labelnew));
  if isfield(input, 'chanunit') && istrue(inverse)
    warning('copying input chanunit to montage');
    [sel1, sel2] = match_str(montage.labelnew, input.label);
    montage.chanunitnew(sel1) = input.chanunit(sel2);
  end
end

if ~isfield(input, 'label') && isfield(input, 'labelnew')
  % the input data structure is also a montage
  inputlabel = input.labelnew;
  if isfield(input, 'chantypenew')
    inputchantype = input.chantypenew;
  else
    inputchantype = repmat({'unknown'}, size(input.labelnew));
  end
  if isfield(input, 'chanunitnew')
    inputchanunit = input.chanunitnew;
  else
    inputchanunit = repmat({'unknown'}, size(input.labelnew));
  end
else
  % the input should describe the channel labels, and optionally the type and unit
  inputlabel = input.label;
  if isfield(input, 'chantype')
    inputchantype = input.chantype;
  else
    inputchantype = repmat({'unknown'}, size(input.label));
  end
  if isfield(input, 'chanunit')
    inputchanunit = input.chanunit;
  else
    inputchanunit = repmat({'unknown'}, size(input.label));
  end
end

% check the consistency of the montage
if ~iscell(montage.labelold) || ~iscell(montage.labelnew)
  error('montage labels need to be specified in cell-arrays');
end

% check the consistency of the montage
if ~all(isfield(montage, {'tra', 'labelold', 'labelnew'}))
  error('the second input argument does not correspond to a montage');
end

% check the consistency of the montage
if size(montage.tra,1)~=length(montage.labelnew)
  error('the number of channels in the montage is inconsistent');
elseif size(montage.tra,2)~=length(montage.labelold)
  error('the number of channels in the montage is inconsistent');
end

% use a default unit transfer from sensors to channels if not otherwise specified
if ~isfield(input, 'tra') && isfield(input, 'label')
  if     isfield(input, 'elecpos') && length(input.label)==size(input.elecpos, 1)
    nchan = length(input.label);
    input.tra = eye(nchan);
  elseif isfield(input, 'coilpos') && length(input.label)==size(input.coilpos, 1)
    nchan = length(input.label);
    input.tra = eye(nchan);
  elseif isfield(input, 'chanpos') && length(input.label)==size(input.chanpos, 1)
    nchan = length(input.label);
    input.tra = eye(nchan);
  end
end

if istrue(inverse)
  % swap the role of the original and new channels
  tmp.labelnew    = montage.labelold;
  tmp.labelold    = montage.labelnew;
  tmp.chantypenew = montage.chantypeold;
  tmp.chantypeold = montage.chantypenew;
  tmp.chanunitnew = montage.chanunitold;
  tmp.chanunitold = montage.chanunitnew;
  % apply the inverse montage, this can be used to undo a previously
  % applied montage
  tmp.tra = full(montage.tra);
  if rank(tmp.tra) < length(tmp.tra)
    warningfun('the linear projection for the montage is not full-rank, the resulting data will have reduced dimensionality');
    tmp.tra = pinv(tmp.tra);
  else
    tmp.tra = inv(tmp.tra);
  end
  montage = tmp;
end

% select and keep the columns that are non-empty, i.e. remove the empty columns
selcol           = find(~all(montage.tra==0, 1));
montage.tra      = montage.tra(:,selcol);
montage.labelold = montage.labelold(selcol);
montage.chantypeold = montage.chantypeold(selcol);
montage.chanunitold = montage.chanunitold(selcol);
clear selcol

% select and remove the columns corresponding to channels that are not present in the
% original data
remove = setdiff(montage.labelold, intersect(montage.labelold, inputlabel));
selcol = match_str(montage.labelold, remove);
% we cannot just remove the colums, all rows that depend on it should also be removed
selrow = false(length(montage.labelnew),1);
for i=1:length(selcol)
  selrow = selrow & (montage.tra(:,selcol(i))~=0);
end
% convert from indices to logical vector
selcol = indx2logical(selcol, length(montage.labelold));
% remove rows and columns
montage.labelold    = montage.labelold(~selcol);
montage.labelnew    = montage.labelnew(~selrow);
montage.chantypeold = montage.chantypeold(~selcol);
montage.chantypenew = montage.chantypenew(~selrow);
montage.chanunitold = montage.chanunitold(~selcol);
montage.chanunitnew = montage.chanunitnew(~selrow);
montage.tra         = montage.tra(~selrow, ~selcol);
clear remove selcol selrow i

% add columns for channels that are present in the input data but not specified in
% the montage, stick to the original order in the data
[dum, ix]   = setdiff(inputlabel, montage.labelold);
addlabel    = inputlabel(sort(ix));
addchantype = inputchantype(sort(ix));
addchanunit = inputchanunit(sort(ix));
m = size(montage.tra,1);
n = size(montage.tra,2);
k = length(addlabel);
% check for NaNs in unused channels; these will be mixed in with the rest
% of the channels and result in NaNs in the output even when multiplied
% with zeros or identity
if k > 0 && isfield(input, 'trial') % check for raw data now only
  cfg = [];
  cfg.channel = addlabel;
  data_unused = ft_selectdata(cfg, input);
  tmp = cat(1, data_unused.trial{:});
  if any(isnan(tmp(:)))
    error('FieldTrip:NaNsinInputData', ['Your input data contains NaNs in channels that are unused '...
      'in the supplied montage. This would result in undesired NaNs in the '...
      'output data. Please remove these channels from the input data (using '...
      'ft_selectdata) before attempting to apply the montage.']);
  end
end
if istrue(keepunused)
  % add the channels that are not rereferenced to the input and output of the
  % montage
  montage.tra((m+(1:k)),(n+(1:k))) = eye(k);
  montage.labelold    = cat(1, montage.labelold(:), addlabel(:));
  montage.labelnew    = cat(1, montage.labelnew(:), addlabel(:));
  montage.chantypeold = cat(1, montage.chantypeold(:), addchantype(:));
  montage.chantypenew = cat(1, montage.chantypenew(:), addchantype(:));
  montage.chanunitold = cat(1, montage.chanunitold(:), addchanunit(:));
  montage.chanunitnew = cat(1, montage.chanunitnew(:), addchanunit(:));
else
  % add the channels that are not rereferenced to the input of the montage only
  montage.tra(:,(n+(1:k))) = zeros(m,k);
  montage.labelold    = cat(1, montage.labelold(:), addlabel(:));
  montage.chantypeold = cat(1, montage.chantypeold(:), addchantype(:));
  montage.chanunitold = cat(1, montage.chanunitold(:), addchanunit(:));
end
clear addlabel addchantype addchanunit m n k

% determine whether all channels are unique
m = size(montage.tra,1);
n = size(montage.tra,2);
if length(unique(montage.labelnew))~=m
  error('not all output channels of the montage are unique');
end
if length(unique(montage.labelold))~=n
  error('not all input channels of the montage are unique');
end

% determine whether all channels that have to be rereferenced are available
if length(intersect(inputlabel, montage.labelold))~=length(montage.labelold)
  error('not all channels that are required in the montage are available in the data');
end

% reorder the columns of the montage matrix
[selinput, selmontage] = match_str(inputlabel, montage.labelold);
montage.tra            = montage.tra(:,selmontage);
montage.labelold       = montage.labelold(selmontage);
montage.chantypeold    = montage.chantypeold(selmontage);
montage.chanunitold    = montage.chanunitold(selmontage);

% ensure that the montage is double precision
montage.tra = double(montage.tra);

% making the tra matrix sparse will speed up subsequent multiplications, but should
% not result in a sparse matrix
% note that this only makes sense for matrices with a lot of zero elements, for dense
% matrices keeping it full will be much quicker
if size(montage.tra,1)>1 && nnz(montage.tra)/numel(montage.tra) < 0.3
  montage.tra = sparse(montage.tra);
else
  montage.tra = full(montage.tra);
end

% update the channel scaling if the input has different units than the montage expects
if isfield(input, 'chanunit') && ~isequal(input.chanunit, montage.chanunitold)
  scale = ft_scalingfactor(input.chanunit, montage.chanunitold);
  montage.tra = montage.tra * diag(scale);
  montage.chanunitold = input.chanunit;
elseif isfield(input, 'chanunitnew') && ~isequal(input.chanunitnew, montage.chanunitold)
  scale = ft_scalingfactor(input.chanunitnew, montage.chanunitold);
  montage.tra = montage.tra * diag(scale);
  montage.chanunitold = input.chanunitnew;
end

if isfield(input, 'chantype') && ~isequal(input.chantype, montage.chantypeold)
  error('inconsistent chantype in data and montage');
elseif isfield(input, 'chantypenew') && ~isequal(input.chantypenew, montage.chantypeold)
  error('inconsistent chantype in data and montage');
end

if isfield(input, 'labelold') && isfield(input, 'labelnew')
  inputtype = 'montage';
elseif isfield(input, 'tra')
  inputtype = 'sens';
elseif isfield(input, 'trial')
  inputtype = 'raw';
elseif isfield(input, 'fourierspctrm')
  inputtype = 'freq';
else
  inputtype = 'unknown';
end

switch inputtype
  case 'montage'
    % apply the montage on top of the other montage
    if isa(input.tra, 'single')
      % sparse matrices and single precision do not match
      input.tra = full(montage.tra) * input.tra;
    else
      input.tra = montage.tra * input.tra;
    end
    input.labelnew    = montage.labelnew;
    input.chantypenew = montage.chantypenew;
    input.chanunitnew = montage.chanunitnew;
    
  case 'sens'
    % apply the montage to an electrode or gradiometer description
    sens = input;
    clear input
    
    % apply the montage to the inputor array
    if isa(sens.tra, 'single')
      % sparse matrices and single precision do not match
      sens.tra = full(montage.tra) * sens.tra;
    else
      sens.tra = montage.tra * sens.tra;
    end
    
    % The montage operates on the coil weights in sens.tra, but the output channels
    % can be different. If possible, we want to keep the original channel positions
    % and orientations.
    [sel1, sel2] = match_str(montage.labelnew, inputlabel);
    keepchans = length(sel1)==length(montage.labelnew);
    
    if isfield(sens, 'chanpos')
      if keepchans
        sens.chanpos = sens.chanpos(sel2,:);
      else
        if ~isfield(sens, 'chanposold')
          % add a chanposold only if it is not there yet
          sens.chanposold = sens.chanpos;
        end
        sens.chanpos = nan(numel(montage.labelnew),3);
      end
    end
    
    if isfield(sens, 'chanori')
      if keepchans
        sens.chanori = sens.chanori(sel2,:);
      else
        if ~isfield(sens, 'chanoriold')
          sens.chanoriold = sens.chanori;
        end
        sens.chanori = nan(numel(montage.labelnew),3);
      end
    end
    
    sens.label    = montage.labelnew;
    sens.chantype = montage.chantypenew;
    sens.chanunit = montage.chanunitnew;
    
    % keep the
    % original label,
    % type and unit
    % for reference
    if ~isfield(sens, 'labelold')
      sens.labelold = inputlabel;
    end
    if ~isfield(sens, 'chantypeold')
      sens.chantypeold = inputchantype;
    end
    if ~isfield(sens, 'chanunitold')
      sens.chanunitold = inputchanunit;
    end
    
    % keep track of the order of the balancing and which one is the current one
    if istrue(inverse)
      if isfield(sens, 'balance')% && isfield(sens.balance, 'previous')
        if isfield(sens.balance, 'previous') && numel(sens.balance.previous)>=1
          sens.balance.current  = sens.balance.previous{1};
          sens.balance.previous = sens.balance.previous(2:end);
        elseif isfield(sens.balance, 'previous')
          sens.balance.current  = 'none';
          sens.balance          = rmfield(sens.balance, 'previous');
        else
          sens.balance.current  = 'none';
        end
      end
    elseif ~istrue(inverse) && ~isempty(bname)
      
      if isfield(sens, 'balance'),
        % check whether a balancing montage with name bname already exist,
        % and if so, how many
        mnt = fieldnames(sens.balance);
        sel = strmatch(bname, mnt);
        if numel(sel)==0,
          % bname can stay the same
        elseif numel(sel)==1
          % the original should be renamed to 'bname1' and the new one should
          % be 'bname2'
          sens.balance.([bname, '1']) = sens.balance.(bname);
          sens.balance                = rmfield(sens.balance, bname);
          if isfield(sens.balance, 'current') && strcmp(sens.balance.current, bname)
            sens.balance.current = [bname, '1'];
          end
          if isfield(sens.balance, 'previous')
            sel2 = strmatch(bname, sens.balance.previous);
            if ~isempty(sel2)
              sens.balance.previous{sel2} = [bname, '1'];
            end
          end
          bname = [bname, '2'];
        else
          bname = [bname, num2str(length(sel)+1)];
        end
      end
      
      if isfield(sens, 'balance') && isfield(sens.balance, 'current')
        if ~isfield(sens.balance, 'previous')
          sens.balance.previous = {};
        end
        sens.balance.previous = [{sens.balance.current} sens.balance.previous];
        sens.balance.current  = bname;
        sens.balance.(bname)  = montage;
      end
    end
    
    % rename the output variable
    input = sens;
    clear sens
    
  case 'raw';
    % apply the montage to the raw data that was preprocessed using fieldtrip
    data = input;
    clear input
    
    Ntrials = numel(data.trial);
    ft_progress('init', feedback, 'processing trials');
    for i=1:Ntrials
      ft_progress(i/Ntrials, 'processing trial %d from %d\n', i, Ntrials);
      if isa(data.trial{i}, 'single')
        % sparse matrices and single
        % precision do not match
        data.trial{i} = full(montage.tra) * data.trial{i};
      else
        data.trial{i} = montage.tra * data.trial{i};
      end
    end
    ft_progress('close');
    
    data.label    = montage.labelnew;
    data.chantype = montage.chantypenew;
    data.chanunit = montage.chanunitnew;
    
    % rename the output variable
    input = data;
    clear data
    
  case 'freq'
    % apply the montage to the spectrally decomposed data
    freq = input;
    clear input
    
    if strcmp(freq.dimord, 'rpttap_chan_freq')
      siz    = size(freq.fourierspctrm);
      nrpt   = siz(1);
      nchan  = siz(2);
      nfreq  = siz(3);
      output = zeros(nrpt, size(montage.tra,1), nfreq);
      for foilop=1:nfreq
        output(:,:,foilop) = freq.fourierspctrm(:,:,foilop) * montage.tra';
      end
      freq.fourierspctrm = output; % replace the original Fourier spectrum
    elseif strcmp(freq.dimord, 'rpttap_chan_freq_time')
      siz    = size(freq.fourierspctrm);
      nrpt   = siz(1);
      nchan  = siz(2);
      nfreq  = siz(3);
      ntime  = siz(4);
      output = zeros(nrpt, size(montage.tra,1), nfreq, ntime);
      for foilop=1:nfreq
        for toilop = 1:ntime
          output(:,:,foilop,toilop) = freq.fourierspctrm(:,:,foilop,toilop) * montage.tra';
        end
      end
      freq.fourierspctrm = output; % replace the original Fourier spectrum
    else
      error('unsupported dimord in frequency data (%s)', freq.dimord);
    end
    
    freq.label    = montage.labelnew;
    freq.chantype = montage.chantypenew;
    freq.chanunit = montage.chanunitnew;
    
    % rename the output variable
    input = freq;
    clear freq
    
  otherwise
    error('unrecognized input');
end % switch inputtype

% only retain the chantype and/or chanunit if they were present in the input
if ~haschantype
  input = removefields(input, {'chantype', 'chantypeold', 'chantypenew'});
end
if ~haschanunit
  input = removefields(input, {'chanunit', 'chanunitold', 'chanunitnew'});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = indx2logical(x, n)
y = false(1,n);
y(x) = true;

function nowarning(varargin)
return

function s = removefields(s, fn)
for i=1:length(fn)
  if isfield(s, fn{i})
    s = rmfield(s, fn{i});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION use "old/new" instead of "org/new"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function montage = fixmontage(montage)
if isfield(montage, 'labelorg')
  montage.labelold = montage.labelorg;
  montage = rmfield(montage, 'labelorg');
end
if isfield(montage, 'chantypeorg')
  montage.chantypeold = montage.chantypeorg;
  montage = rmfield(montage, 'chantypeorg');
end
if isfield(montage, 'chanunitorg')
  montage.chanunitold = montage.chanunitorg;
  montage = rmfield(montage, 'chanunitorg');
end


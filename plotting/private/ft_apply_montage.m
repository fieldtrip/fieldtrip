function [input] = ft_apply_montage(input, montage, varargin)

% FT_APPLY_MONTAGE changes the montage of an electrode or gradiometer array. A
% montage can be used for EEG rereferencing, MEG synthetic gradients, MEG
% planar gradients or unmixing using ICA. This function applies the montage
% to the inputor array. The inputor array can subsequently be used for
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
%   montage.labelnew = Mx1 cell-array
%   montage.labelorg = Nx1 cell-array
%
% As an example, a bipolar montage could look like this
%   bipolar.labelorg  = {'1', '2', '3', '4'}
%   bipolar.labelnew  = {'1-2', '2-3', '3-4'}
%   bipolar.tra       = [
%     +1 -1  0  0
%      0 +1 -1  0
%      0  0 +1 -1
%   ];
%
% Additional options should be specified in key-value pairs and can be
%   'keepunused'    string, 'yes' or 'no' (default = 'no')
%   'inverse'       string, 'yes' or 'no' (default = 'no')
%
% If the first input is a montage, then the second input montage will be
% applied to the first. In effect, the output montage will first do
% montage1, then montage2.
%
% See also FT_READ_SENS, FT_TRANSFORM_SENS

% Copyright (C) 2008-2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

% get optional input arguments
keepunused  = ft_getopt(varargin, 'keepunused',  'no');
inverse     = ft_getopt(varargin, 'inverse',     'no');
feedback    = ft_getopt(varargin, 'feedback',    'text');
showwarning = ft_getopt(varargin, 'warning',     'yes');
bname       = ft_getopt(varargin, 'balancename', '');

if ~isfield(input, 'label') && isfield(input, 'labelnew')
  % the input data structure is also a montage
  inputlabel = input.labelnew;
else
  % the input should describe the channel labels
  inputlabel = input.label;
end

if strcmp(showwarning, 'yes')
  warningfun = @warning;
else
  warningfun = @nowarning;
end

% check the consistency of the input inputor array or data
if ~all(isfield(montage, {'tra', 'labelorg', 'labelnew'}))
  error('the second input argument does not correspond to a montage');
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
    warningfun('the linear projection for the montage is not full-rank, the resulting data will have reduced dimensionality');
    tmp.tra = pinv(tmp.tra);
  else
    tmp.tra = inv(tmp.tra);
  end
  montage = tmp;
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

% select and keep the columns that are non-empty, i.e. remove the empty columns
selcol           = find(~all(montage.tra==0, 1));
montage.tra      = montage.tra(:,selcol);
montage.labelorg = montage.labelorg(selcol);
clear selcol

% select and remove the columns corresponding to channels that are not present in the original data
remove = setdiff(montage.labelorg, intersect(montage.labelorg, inputlabel));
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
[add, ix] = setdiff(inputlabel, montage.labelorg);
add = inputlabel(sort(ix));
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
if length(intersect(inputlabel, montage.labelorg))~=length(montage.labelorg)
  error('not all channels that are required in the montage are available in the data');
end

% reorder the columns of the montage matrix
[selinput, selmontage]  = match_str(inputlabel, montage.labelorg);
montage.tra             = double(montage.tra(:,selmontage));
montage.labelorg        = montage.labelorg(selmontage);

% making the tra matrix sparse will speed up subsequent multiplications
% but should not result in a sparse matrix
% note that this only makes sense for matrices with a lot of zero elements,
% for dense matrices keeping it full will be much quicker
if size(montage.tra,1)>1 && nnz(montage.tra)/numel(montage.tra) < 0.3
  montage.tra = sparse(montage.tra);
end

inputtype = 'unknown';
if isfield(input, 'labelorg') && isfield(input, 'labelnew')
  inputtype = 'montage';
elseif isfield(input, 'tra')
  inputtype = 'sens';
elseif isfield(input, 'trial')
  inputtype = 'raw';
elseif isfield(input, 'fourierspctrm')
  inputtype = 'freq';
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
    input.labelnew = montage.labelnew;
    
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
    
    % The montage operates on the coil weights in sens.tra, but the output
    % channels can be different. If possible, we want to keep the original
    % channel positions and orientations.
    [sel1, sel2] = match_str(montage.labelnew, inputlabel);
    keepchans = length(sel1)==length(montage.labelnew);
    
    if isfield(sens, 'chanpos')
      if keepchans
        sens.chanpos = sens.chanpos(sel2,:);
      else
        sens.chanpos = nan(numel(montage.labelnew),3);
      end
    end
    
    if isfield(sens, 'chanori')
      if keepchans
        sens.chanori = sens.chanori(sel2,:);
      else
        sens.chanori = nan(numel(montage.labelnew),3);
      end
    end
    
    if isfield(sens, 'chantype')
      if keepchans
        sens.chantype = sens.chantype(sel2,:);
      else
        % FIXME don't know what to do here  
        sens = rmfield(sens, 'chantype');
      end
    end
    
    if isfield(sens, 'chanunit')
      if keepchans
        sens.chanunit = sens.chanunit(sel2,:);
      else
        % FIXME don't know what to do here  
        sens = rmfield(sens, 'chanunit');
      end
    end
    
    sens.label = montage.labelnew;
    
    % keep track of the order of the balancing and which one is the current one
    if strcmp(inverse, 'yes')
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
    elseif ~strcmp(inverse, 'yes') && ~isempty(bname)
      
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
        % sparse matrices and single precision do not match
        data.trial{i} = full(montage.tra) * data.trial{i};
      else
        data.trial{i} = montage.tra * data.trial{i};
      end
    end
    ft_progress('close');
    
    data.label = montage.labelnew;
    
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
    else
      error('unsupported dimord in frequency data (%s)', freq.dimord);
    end
    
    % replace the Fourier spectrum
    freq.fourierspctrm = output;
    freq.label = montage.labelnew;
    
    % rename the output variable
    input = freq;
    clear freq
    
  otherwise
    error('unrecognized input');
end % switch inputtype

% check whether the input contains chantype and/or chanunit and remove these
% as they may have been invalidated by the transform (e.g. with megplanar)
[sel1, sel2] = match_str(montage.labelnew, inputlabel);
keepchans = (length(sel1)==length(montage.labelnew));

if isfield(input, 'chantype')
  if keepchans
    % reorder them according to the montage
    sens.chantype = input.chantype(sel2,:);
  else
    input = rmfield(input, 'chantype');
  end
end

if isfield(input, 'chanunit')
  if keepchans
    % reorder them according to the montage
    sens.chanunit = input.chanunit(sel2,:);
  else
    input = rmfield(input, 'chanunit');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = indx2logical(x, n)
y = false(1,n);
y(x) = true;

function nowarning(varargin)
return


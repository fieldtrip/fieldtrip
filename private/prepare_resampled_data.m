function [cfg, varargout] = prepare_resampled_data(cfg, varargin)

% PREPARE_RESAMPLED_DATA performs resampling of the input data for
% multiple variables in a single or multiple conditions. The resampling
% will be performed along the first dimension of every input variable. This
% function is intended to be used as subfunction for various algorithms
% implemented in FieldTrip.
% 
% Supported resampling strategies are
%   jackknife for one condition
%   bootstrap for one condition
%   permutation for two conditions
%   resampling for two or more conditions
% You can also specify that you do not want any resampling, in which case
% only the average over the original data will be computed.
% 
% Use as
%   [cfg, varargout] = prepare_resampled_data(cfg, varargin)
% where the configuration can contain
%   cfg.jackknife        = 'yes' or 'no'
%   cfg.bootstrap        = 'yes' or 'no'
%   cfg.pseudovalue      = 'yes' or 'no'
%   cfg.randomization    = 'yes' or 'no'
%   cfg.permutation      = 'yes' or 'no'
%   cfg.numbootstrap     = number
%   cfg.numrandomization = number
%   cfg.numpermutation   = number, or 'all'
% and the input and output data is orgainzed according to the examples below.
% 
% for N data objects in one condition
%   [cfg, r1, r2 ... rN] = prepare_resampled_data(cfg, o1, o2 ... oN) 
%
% for N data objects in two conditions
%   [cfg, r11 ... r1N, r21 ... rN] = prepare_resampled_data(cfg, o11 ... o1N, o21 ... o2N)
%
% for multiple data objects in three conditions
%   [cfg, r11..., r21 ..., r31 ...] = prepare_resampled_data(cfg, o11 ..., o21 ..., o31 ...);

% For definitions (not used here!) and for other resampling approaches see
% http://www.nist.gov/dads/HTML/permutation.html

% The definitions used here are
% PERMUTATION a pair-wise rearrangement of elements over two categories
% RANDOMIZATION a rearrangement of elements over two categories, where
%   the number of elements in each category remains the same

% Copyright (C) 2004, Robert Oostenveld
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

% for backward compatibility with misspelled configuration option
ft_checkconfig(cfg, 'renamed', {'jacknife', 'jackknife'});

% set the defaults
if ~isfield(cfg, 'jackknife'),        cfg.jackknife        = 'no';  end
if ~isfield(cfg, 'bootstrap'),        cfg.bootstrap        = 'no';  end
if ~isfield(cfg, 'pseudovalue'),      cfg.pseudovalue      = 'no';  end
if ~isfield(cfg, 'randomization'),    cfg.randomization    = 'no';  end
if ~isfield(cfg, 'permutation'),      cfg.permutation      = 'no';  end
if ~isfield(cfg, 'numbootstrap'),     cfg.numbootstrap     = 100;   end
if ~isfield(cfg, 'numrandomization'), cfg.numrandomization = 100;   end
if ~isfield(cfg, 'numpermutation'),   cfg.numpermutation   = 'all'; end

% remove configuration fields that are not needed
if strcmp(cfg.bootstrap, 'no'), cfg = rmfield(cfg, 'numbootstrap'); end
if strcmp(cfg.randomization, 'no'), cfg = rmfield(cfg, 'numrandomization'); end
if strcmp(cfg.permutation, 'no'), cfg = rmfield(cfg, 'numpermutation'); end

% use some boolean flags for convenience
resampleaverage       = 0;
resamplejackknife      = strcmp(cfg.jackknife, 'yes');
resamplebootstrap     = strcmp(cfg.bootstrap, 'yes');
resamplepseudovalue   = strcmp(cfg.pseudovalue, 'yes');
resamplerandomization = strcmp(cfg.randomization, 'yes');
resamplepermutation   = strcmp(cfg.permutation, 'yes');

tmp = resamplejackknife+resamplebootstrap+resamplepseudovalue+resamplerandomization+resamplepermutation;
if tmp>1
  error('only one resampling strategy should be specified');
elseif tmp<1
  % default is to average over trials
  resampleaverage = 1;
end

if ~isfield(cfg, 'numcondition')
  % set the default number of conditions for each method
  if resampleaverage || resamplejackknife || resamplebootstrap || resamplepseudovalue
    % FIXME for averaging the number of conditions can also be larger than 1
    cfg.numcondition = 1;
  elseif resamplerandomization || resamplepermutation
    % FIXME for randomization the number of conditions can already be larger than 2
    % and for permutation it in principle can also be larger than two
    cfg.numcondition = 2;
  end
end

% distribute the input dataobjects over the conditions
Ncondition = cfg.numcondition;
Ndata = (nargin-1)/Ncondition; % the first argument is the cfg
if Ndata~=round(Ndata)
  % this should not be fractional
  error('incorrect number of input arguments');
else
  for c=1:Ncondition
    for d=1:Ndata
      datain{c,d} = varargin{(c-1)*Ndata + d};
    end
  end
end

fprintf('%d conditions, each with %d data objects\n', Ncondition, Ndata);

for c=1:Ncondition
  for d=1:Ndata
    % the first dimension encodes replications
    siz = size(datain{c,d});
    Nreplication(c,d) = siz(1);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if resampleaverage
  % the number of replications of all dataobjects should be the same 
  for c=1:Ncondition
    fprintf('averaging %d replications for one condition\n', max(Nreplication(c,:)));
    for d=1:Ndata
      siz = size(datain{c,d});
      if siz(1)==max(Nreplication(c,:))
        dataout{c,d} = reshape(mean(datain{c,d}, 1), [siz(2:end) 1]);
      elseif siz(1)==1
        dataout{c,d} = datain{c,d};
      else
        error('inconsistent number of replications in data');
      end
    end 
  end 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif resamplejackknife
  if Ncondition>1
    error('jackknife requires exactly one condition');
  end
  % the number of replications of all dataobjects should be the same 
  Nreplication = max(Nreplication(:));
  fprintf('constructing %d jackknife replications\n', Nreplication);
  for d=1:Ndata
    siz = size(datain{1,d});
    dataout{1,d} = zeros(siz);     % output dimensions are the same as input
    if siz(1)==Nreplication
      % compute the sum over all replications
      s = sum(datain{1,d}, 1);
      for r=1:Nreplication
        % compute the average, leaving this replication out
        dataout{1,d}(r,:) = (s(:)' - datain{1,d}(r,:))./(Nreplication-1);
      end
    elseif siz(1)==1
      dataout{1,d} = datain{1,d};
    else
      error('inconsistent number of replications in data');
    end
  end 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif resamplebootstrap
  if Ncondition>1
    error('bootstrap requires exactly one condition');
  end
  % the number of replications of all dataobjects should be the same 
  Nreplication = max(Nreplication(:));
  fprintf('constructing %d bootstrap replications\n', cfg.numbootstrap);
  % create the random bootstrap replications
  bootstrap = ceil(Nreplication*rand(cfg.numbootstrap, Nreplication));
  for d=1:Ndata
    siz = size(datain{1,d});
    dataout{1,d} = zeros([cfg.numbootstrap siz(2:end)]);
    if siz(1)==Nreplication
      for r=1:cfg.numbootstrap
        % compute the average for this bootstrap replication
        dataout{1,d}(r,:) = mean(datain{1,d}(bootstrap(r,:),:), 1);
      end
    elseif siz(1)==1
      dataout{1,d} = datain{1,d};
    else
      error('inconsistent number of replications in data');
    end
  end 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif resamplepseudovalue
  if Ncondition>1
    error('pseudovalue requires exactly one condition');
  end
  % the number of replications of all dataobjects should be the same 
  Nreplication = max(Nreplication(:));
  fprintf('constructing %d pseudovalue replications\n', Nreplication);
  for d=1:Ndata
    siz = size(datain{1,d});
    siz(1) = siz(1) + 1;
    dataout{1,d} = zeros(siz);     % output dimensions are the same as input, except for one extra replication
    if siz(1)==Nreplication+1
      % compute the sum over all replications
      s = sum(datain{1,d}, 1);
      % compute the average over all replications
      dataout{1,d}(1,:) = s(:)'/Nreplication;
      for r=1:Nreplication
        % compute the average, leaving this replication out
        dataout{1,d}(r+1,:) = (s(:)' - datain{1,d}(r,:))./(Nreplication-1);
      end
    elseif siz(1)==1
      error('multiple replications are required in data');
    else
      error('inconsistent number of replications in data');
    end
  end 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif resamplerandomization
  if Ncondition<2
    error('randomization requires at least two conditions');
  end
  
  % the number of replications should be the same for the different dataobjects 
  % within one condition, but can be different over conditions
  tmp = [];
  for c=1:Ncondition
    tmp(c) = max(Nreplication(c,:));
  end
  Nreplication = tmp(:);
  
  fprintf('creating %d random resamplings\n', cfg.numrandomization);
  for r=1:cfg.numrandomization
    [dum, randomlist] = sort(rand(1,sum(Nreplication)));
    resample(r,:) = randomlist;
  end
  
  for d=1:Ndata
    siz = [];
    for c=1:Ncondition
      siz(c,:) = size(datain{c,d});
    end
    for c=2:Ncondition
      if ~all(siz(c, 2:end)==siz(1, 2:end))
        error('the dimensions should be the same for every replication in every condition')
      end
    end
    % concatenate the data in all conditions to facilitate the random resamplings
    dataconcat = cat(1, datain{:,d});
    % allocate memory to hold the output
    for c=1:Ncondition
      dataout{c,d} = zeros([cfg.numrandomization siz(1, 2:end)]);
    end
    if all(siz(:,1)==Nreplication)
      for r=1:cfg.numrandomization
        for c=1:Ncondition
          selbeg = sum(Nreplication(1:(c-1))) + 1;
          selend = selbeg + Nreplication(c) - 1;
          % compute the average for this random replication
          dataout{c,d}(r,:) = mean(dataconcat(resample(r,selbeg:selend),:), 1);
        end
      end
    elseif all(siz(:,1)==1)
      for c=1:Ncondition
        dataout{c,d} = datain{c,d};
      end
    else
      error('inconsistent number of replications in data');
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif resamplepermutation
  if Ncondition<2
    error('permutation requires at least two conditions');
  end
  
  % the number of replications should be the same for all dataobjects in both conditions
  Nreplication = max(Nreplication(:));
  
  if Ncondition<=2
    % first create all possible permutations and then make a random
    % selection out of it (or select all of them in case of 'all')
    
    if ischar(cfg.numpermutation) && strcmp(cfg.numpermutation, 'all')
      % create all possible permutations
      cfg.numpermutation = 2^Nreplication;
      fprintf('creating all (%d) permutations\n', cfg.numpermutation);
      permute = char(zeros(cfg.numpermutation, Nreplication));
      for i=1:cfg.numpermutation
        permute(i,:) = dec2bin((i-1), Nreplication);
      end
      % convert char array into matrix with booleans
      permute = (permute=='1');
    else
      % only create a selection of all possible permutations to save memory
      fprintf('creating %d random permutations from total %d\n', cfg.numpermutation, 2^Nreplication);
      % The idea is that every binary number can represent a permutation,
      % therefore we have to make a list of random numbers between one and
      % the maximum number of permutations. Normally I would make a list with
      % all numbers and then make a random selection out of that list, but the
      % maximum number of permutations can be extremely large. Therefore I make
      % a list of unique random numbers, and keep extending that list untill I
      % have enough unique numbers.
      if cfg.numpermutation>2^Nreplication
        error(sprintf('it is not possible to make more random permutations than %d', 2^Nreplication));
      end
      number = [];
      while length(number)<cfg.numpermutation
        % iterate to obtain enough unique random numbers
        number = [number round(rand(1, cfg.numpermutation)*(2^Nreplication-1) + 1)];
        number = unique(number);
      end
      if length(number)>cfg.numpermutation
        % we got too many, remove some of the permutations
        randomlist = randperm(length(number));
        number = number(randomlist(1:cfg.numpermutation));
      end
      % convert the list with random numbers into their binary representations
      permute = char(zeros(cfg.numpermutation, Nreplication));
      for i=1:cfg.numpermutation
        permute(i,:) = dec2bin((number(i)-1), Nreplication);
      end
      % convert char array into matrix with booleans
      permute = (permute=='1');
    end
    
    % convert the binary/boolean matrix that only indicates for each trial whether it should be permuted
    % into a selection array that can be used as indices into the concatenated data
    tmp = [];
    for r=1:cfg.numpermutation
      sel1 = (1:Nreplication) + (permute(r,:)==0)*Nreplication; % select the trials assigned to condition 1
      sel2 = (1:Nreplication) + (permute(r,:)==1)*Nreplication; % select the trials assigned to condition 2
      tmp(r,:) = [sel1 sel2];
    end 
    permute = tmp;
    
  else
    % Alternative implementation for more than two conditions, where it is not
    % possible to generate all possible permutations.  The idea is that every
    % permutation of one trial can be written as a number between 1 and the
    % maximum number of possible permutations.  Hence create a matrix with all
    % permutations, and for every trial select one of these permutations. Then
    % combine this permutation over all trials, and repeat untill theer are
    % enough resamplings of all trials.
    if ischar(cfg.numpermutation) && strcmp(cfg.numpermutation, 'all')
      error('it is not possible to generate all possible permutations for more than two conditions');
    end
    allperms = sortrows(perms(1:Ncondition));
    numperms = prod(1:Ncondition);      % equal to size(allperms,1)
    permute = zeros(cfg.numpermutation, Ncondition*Nreplication);
    for r=1:cfg.numpermutation
      % select one of the possible permutations for every trial
      sel = floor(rand(1,Nreplication)*numperms) + 1;   
      tmp = allperms(sel,:);
      % convert to a list that can be used as indices into the concatenated data
      tmp = (tmp-1)*Nreplication + repmat((1:Nreplication)', [1 Ncondition]);
      permute(r,:) = tmp(:)';
    end
  end
  
  for d=1:Ndata
    % determine the size of this data object 
    siz = size(datain{1,d});
    for c=2:Ncondition
      if ~all(siz==size(datain{c,d}))
        error('the dimensions should be the same for every condition')
      end
    end
    % concatenate the data in all conditions to facilitate the random resamplings
    dataconcat = cat(1, datain{:,d});
    % allocate memory to hold the output
    for c=1:Ncondition
      dataout{c,d} = zeros([cfg.numpermutation siz(2:end)]);
    end
    if siz(1)==Nreplication
      for r=1:cfg.numpermutation
        for c=1:Ncondition
          sel = permute(r, ((c-1)*Nreplication+1):(c*Nreplication));
          dataout{c,d}(r,:) = mean(dataconcat(sel,:), 1);
        end
      end
    elseif siz(1)==1
      for c=1:Ncondition
        dataout{c,d} = datain{c,d};
      end
    else
      error('inconsistent number of replications in data');
    end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  error('unknown resampling strategy')
end

% reassign the resampled dataobjects over the output arguments
for c=1:Ncondition
  for d=1:Ndata
    varargout{(c-1)*Ndata + d} = dataout{c,d};
  end
end


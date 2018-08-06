function output = ft_connectivity_mutualinformation(input, varargin)

% FT_CONNECTIVITY_MUTUALINFORMATION computes mutual information using the information
% breakdown toolbox (ibtb), as described in Magri et al., BMC Neuroscience 2009,
% 1471-2202.
%
% Use as
%   mi = ft_connectivity_mutualinformation(data, ...)
%
% The input data should be a Nchan x Nobservations matrix.
%
% Additional optional input arguments come as key-value pairs:
%   histmethod = The way that histograms are generated from the data. Possible values
%                are 'eqpop' (default), 'eqspace', 'ceqspace', 'gseqspace'.
%                See the help of the 'binr' function in the ibtb toolbox for more information.
%   numbin     = scalar value. The number of bins used to create the histograms needed for
%                the entropy computations
%   opts       = structure that is passed on to the 'information' function in the ibtb
%                toolbox. See the help of that function for more information.
%   refindx    = scalar value or 'all'. The channel that is used as 'reference channel'.
%   refdata    = 1xNobservations vector, as an alternative to the refindx. Refdata takes precedence over refindx
%
% The output contains the estimated mutual information between all channels and the reference channel(s).
%
% See also FT_CONNECTIVITYANALYSIS

% Copyright (C) 2016 Donders Institute, Jan-Mathijs Schoffelen
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


% check whether the required toolbox is available
ft_hastoolbox('ibtb', 1);

% set some options
histmethod  = ft_getopt(varargin, 'histmethod', 'eqpop');
numbin      = ft_getopt(varargin, 'numbin',     10);
refindx     = ft_getopt(varargin, 'refindx',    'all');
refdata     = ft_getopt(varargin, 'refdata',    []);
lags        = ft_getopt(varargin, 'lags',       0); % shift of data w.r.t. reference, in samples

% set some additional options that pertain to the algorithmic details of the mutual information computation
opts        = ft_getopt(varargin, 'opts',       []);
opts.nt     = ft_getopt(opts,     'nt',         []);
opts.method = ft_getopt(opts,     'method',     'dr');
opts.bias   = ft_getopt(opts,     'bias',       'pt');

if ischar(refindx) && strcmp(refindx, 'all')
  refindx = (1:size(input,1))';
end

if numel(lags)>1 || lags~=0
  if numel(refindx)>1, ft_error('with multiple lags, or with a lag~=0 only a single refindx is allowed'); end
  refdata = input(refindx,:);
  n       = size(refdata,2);
  
  output = zeros(size(input,1), numel(lags));
  for k = 1:numel(lags)
    fprintf('computing mutualinformation for time lag in samples %d\n', lags(k));
    
    beg1 = max(0, lags(k))  + 1;
    beg2 = max(0, -lags(k)) + 1;
    n1   = n-abs(lags(k));
    
    end1 = beg1+n1-1;
    end2 = beg2+n1-1;
    
    tmpdata = nan(size(refdata));
    tmpdata(beg1:end1) = refdata(beg2:end2);
    
    tmp = ft_connectivity_mutualinformation(input, 'refdata', tmpdata, 'opts', opts, 'histmethod', histmethod, 'numbin', numbin);
    output(:,k) = tmp(1:end-1);
  end
  return
end

if ~isempty(refdata)
  input   = cat(1, input, refdata);
  refindx = size(input,1);
end

% check validity of refindx
if length(refindx)~=numel(refindx)
  % could be channelcmb indexing
  ft_error('channelcmb indexing is not supported');
end

% get rid of nans in the input
notsel = sum(~isfinite(input))>0;
input  = input(:,~notsel);

[nchan, nsmp] = size(input);
output        = zeros(nchan, numel(refindx))+nan;

for k = 1:numel(refindx)
  signal1 = input(refindx(k),:);
  
  % discretize signal1
  signal1 = binr(signal1, nsmp, numbin, histmethod);
  
  for m = setdiff(1:size(input,1),refindx(k))
    signal2 = input(m,:);
    
    % represent signal2 in bins according to signal1's discretization
    R = zeros(1,3,numbin);
    for j = 1:numbin
      nr         = signal1==j-1;
      opts.nt(j) = sum(nr);
      R(1, 1:opts.nt(j),j) = signal2(nr);
    end
    
    % discretize signal2 and compute mi
    R2 = binr(R, opts.nt', numbin, histmethod);
    output(m,k) = information(R2, opts, 'I'); % this computes mutual information
    
  end
end

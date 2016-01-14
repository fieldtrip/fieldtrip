function output = ft_connectivity_mutualinformation(input, varargin)

% FT_CONNECTIVITY_MUTUALINFORMATION computes mutual information using 
% the information breakdown toolbox (ibtb), as described in Magri et al.,
% BMC Neuroscience 2009, 1471-2202. The function is a helper function for
% FT_CONNECTIVITYANALYSIS. As a standalone function, it could be used as
% follows:
%
% mi = ft_connectivity_mutualinformation(data, varargin)
%
% Additional input arguments come as key-value pairs:
%   histmethod = The way that histograms are generated from the data. Possible values
%                are 'eqpop' (default), 'eqspace', 'ceqspace', 'gseqspace'.
%                See the help of the 'binr' function in the ibtb toolbox for more information.
%   numbin     = scalar value. The number of bins used to create the histograms needed for 
%                the entropy computations
%   refindx    = scalar value or 'all'. The channel that is used as 'reference channel'.
%   opts       = structure that is passed on to the 'information' function in the ibtb
%                toolbox. See the help of that function for more information.
%
% The output contains the estimated mutual information between all channels and the reference channel(s).

% Copyright (C) 2016 Donders Institute, Jan-Mathijs Schoffelen
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
% $Id: ft_connectivity_corr.m 9965 2014-11-13 16:56:13Z roboos $

 
% check whether the required toolbox is available
ft_hastoolbox('ibtb', 1);

% set some options
histmethod = ft_getopt(varargin, 'histmethod', 'eqpop');
numbin     = ft_getopt(varargin, 'numbin',     10);
refindx    = ft_getopt(varargin, 'refindx',    'all');

if ischar(refindx) && strcmp(refindx, 'all')
  refindx = (1:size(input,1))';
end

% check validity of refindx
if length(refindx)~=numel(refindx)
  % could be channelcmb indexing
  error('channelcmb indexing is not, yet supported');
end

% set some additional options that pertain to the algorithmic details of the
% mutual information computation
opts        = ft_getopt(varargin, 'opts', []);
opts.nt     = ft_getopt(opts, 'nt', []);
opts.method = ft_getopt(opts, 'method', 'dr');
opts.bias   = ft_getopt(opts, 'bias',   'pt');

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

function [c, v, n] = ft_connectivity_psi(input, varargin)

% FT_CONNECTIVITY_PSI computes the phase slope index from a data-matrix
% containing the cross-spectral density. It implements the method described
% in: Nolte et al., Robustly estimating the flow direction of information
% in complex physical systems. Physical Review Letters, 2008; 100; 234101.
%
% Use as
%   [c, v, n] = ft_connectivity_psi(input, ...)
%
% The input data input should be organized as
%   Repetitions x Channel x Channel (x Frequency) (x Time)
% or
%   Repetitions x Channelcombination (x Frequency) (x Time)
%
% The first dimension should be singleton if the input already contains an
% average.
%
% Additional optional input arguments come as key-value pairs:
%   nbin			=	scalar, half-bandwidth parameter: the number of frequency bins
%								across which to integrate
%   hasjack		= 0 or 1, specifying whether the repetitions represent
%               leave-one-out samples (allowing for a variance estimate)
%   feedback	= 'none', 'text', 'textbar' type of feedback showing progress of
%               computation
%   dimord		= string, specifying how the input matrix should be interpreted
%   powindx   =
%   normalize =
%
% The output p contains the phase slope index, v is a variance estimate
% which only can be computed if the data contains leave-one-out samples,
% and n is the number of repetitions in the input data. If the phase slope
% index is positive, then the first chan (1st dim) becomes more lagged (or
% less leading) with higher frequency, indicating that it is causally
% driven by the second channel (2nd dim)
%
% See also FT_CONNECTIVITYANALYSIS

% Copyright (C) 2009-2010 Donders Institute, Jan-Mathijs Schoffelen
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

% FIXME: interpretation of the slope

hasjack   = ft_getopt(varargin, 'hasjack', 0);
feedback  = ft_getopt(varargin, 'feedback', 'none');
dimord    = ft_getopt(varargin, 'dimord');
powindx   = ft_getopt(varargin, 'powindx');
normalize = ft_getopt(varargin, 'normalize', 'no');
nbin      = ft_getopt(varargin, 'nbin');

if isempty(dimord)
  ft_error('input parameters should contain a dimord');
end

if (length(strfind(dimord, 'chan'))~=2 || contains(dimord, 'pos')>0) && ~isempty(powindx)
  %crossterms are not described with chan_chan_therest, but are linearly indexed
  
  siz = size(input);
  
  outsum = zeros(siz(2:end));
  outssq = zeros(siz(2:end));
  pvec   = [2 setdiff(1:numel(siz),2)];
  
  ft_progress('init', feedback, 'computing metric...');
  %first compute coherency and then phaseslopeindex
  for j = 1:siz(1)
    ft_progress(j/siz(1), 'computing metric for replicate %d from %d\n', j, siz(1));
    c      = reshape(input(j,:,:,:,:), siz(2:end));
    p1     = abs(reshape(input(j,powindx(:,1),:,:,:), siz(2:end)));
    p2     = abs(reshape(input(j,powindx(:,2),:,:,:), siz(2:end)));
    
    p      = ipermute(phaseslope(permute(c./sqrt(p1.*p2), pvec), nbin, normalize), pvec);
    
    outsum = outsum + p;
    outssq = outssq + p.^2;
  end
  ft_progress('close');
  
elseif length(strfind(dimord, 'chan'))==2 || length(strfind(dimord, 'pos'))==2
  %crossterms are described by chan_chan_therest
  
  siz = size(input);
  
  outsum = zeros(siz(2:end));
  outssq = zeros(siz(2:end));
  pvec   = [3 setdiff(1:numel(siz),3)];
  
  ft_progress('init', feedback, 'computing metric...');
  for j = 1:siz(1)
    ft_progress(j/siz(1), 'computing metric for replicate %d from %d\n', j, siz(1));
    p1  = zeros([siz(2) 1 siz(4:end)]);
    p2  = zeros([1 siz(3) siz(4:end)]);
    for k = 1:siz(2)
      p1(k,1,:,:,:,:) = input(j,k,k,:,:,:,:);
      p2(1,k,:,:,:,:) = input(j,k,k,:,:,:,:);
    end
    c      = reshape(input(j,:,:,:,:,:,:), siz(2:end));
    p1     = p1(:,ones(1,siz(3)),:,:,:,:);
    p2     = p2(ones(1,siz(2)),:,:,:,:,:);
    p      = ipermute(phaseslope(permute(c./sqrt(p1.*p2), pvec), nbin, normalize), pvec);
    p(isnan(p)) = 0;
    outsum = outsum + p;
    outssq = outssq + p.^2;
  end
  ft_progress('close');
  
end

n = siz(1);
c = outsum./n;

if n>1
  n = shiftdim(sum(~isnan(input),1),1);
  if hasjack
    bias = (n-1).^2;
  else
    bias = 1;
  end
  v = bias.*(outssq - (outsum.^2)./n)./(n - 1);
else
  v = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y] = phaseslope(x, n, norm)

m   = size(x, 1); %total number of frequency bins
y   = zeros(size(x));
x(1:end-1,:,:,:,:) = conj(x(1:end-1,:,:,:,:)).*x(2:end,:,:,:,:);

if strcmp(norm, 'yes')
  coh = zeros(size(x));
  coh(1:end-1,:,:,:,:) = (abs(x(1:end-1,:,:,:,:)) .* abs(x(2:end,:,:,:,:))) + 1;
  %FIXME why the +1? get the coherence
  for k = 1:m
    begindx = max(1,k-n);
    endindx = min(m,k+n);
    y(k,:,:,:,:) = imag(nansum(x(begindx:endindx,:,:,:,:)./coh(begindx:endindx,:,:,:,:),1));
  end
else
  for k = 1:m
    begindx = max(1,k-n);
    endindx = min(m,k+n);
    y(k,:,:,:,:) = imag(nansum(x(begindx:endindx,:,:,:,:),1));
  end
end

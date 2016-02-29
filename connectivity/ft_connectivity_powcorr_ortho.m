function [c] = ft_connectivity_powcorr_ortho(mom, varargin)

% FT_CONNECTIVITY_POWCORR_ORTHO computes power correlation after removing
% the zero-lag contribution on a trial-by-trial basis, according to Hipp's
% Nature Neuroscience paper. 
%
% Use as
%   c = ft_connectivity_powcorr(mom)
%   c = ft_connectivity_powcorr(mom, 'refindx', refindx)
%
% Where mom is a NchanxNrpt matrix containing the complex-valued amplitude
% and phase information at a given frequency, and the optional key refindx
% specifies the index/indices of the channels that serve as a reference 
% channel. (Default is 'all').
%
% The output c is a NchanxNrefchan matrix that contain the power correlation
% for all channels orthogonalised relative to the reference channel in the first
% Nrefchan columns, and the power correlation for the reference channels 
% orthogonalised relative to the channels in the second Nrefchan columns.

% Copyright (C) 2012 Jan-Mathijs Schoffelen
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

refindx = ft_getopt(varargin, 'refindx', 'all');
tapvec  = ft_getopt(varargin, 'tapvec',  ones(1,size(mom,2)));

if strcmp(refindx, 'all')
  refindx = 1:size(mom,1);
end

cmomnorm = conj(mom./abs(mom)); % only need to do conj() once

n        = size(mom,1);
ntap     = tapvec(1);
if ~all(tapvec==ntap)
  error('unequal number of tapers per observation is not yet supported');
end
% FIXME think about multiple tapers per trial
%if ntap>1
%  error('more than one taper per observation is not yet supported');
%end

% create a sparse matrix tra, that can be used as a right multiplying
% matrix to combine across tapers

ix = zeros(sum(tapvec),1);
jx = ix;
sx = ix;
for k = 1:numel(tapvec)
  indx = (k-1)*ntap+(1:ntap);
  ix(indx) = indx;
  jx(indx) = k;
  sx(indx) = 1./ntap;
end
tra = sparse(ix,jx,sx,sum(tapvec),numel(tapvec));

%tra  = zeros(size(mom,2), numel(tapvec));
%for k = 1:numel(tapvec)
%  tra((k-1)*ntap+(1:ntap), k) = 1./ntap;
%end

powmom = (abs(mom).^2)*tra; % need only once
powmom = standardise(log10(powmom), 2);

c = zeros(n, numel(refindx));%;*2);
N = ones(n,1);
%warning off;
for k = 1:numel(refindx)      
  indx     = refindx(k);
  ref      = mom(indx,:);
  crefnorm = conj(ref./abs(ref));

  % FIXME the following is probably not correct for ntap>1
  pow2 = (abs(imag(ref(N,:).*cmomnorm)).^2)*tra;
  pow2 = standardise(log10(pow2), 2);
  c1   = mean(powmom.*pow2, 2);
  pow1 = (abs(imag(mom.*crefnorm(N,:))).^2)*tra;
  pow1 = standardise(log10(pow1), 2);
  
  pow2 = (abs(ref).^2)*tra;
  pow2 = standardise(log10(pow2), 2);
  pow2 = repmat(pow2, [n 1]);
  c2   = mean(pow1.*pow2, 2);

  c(:,k) = (c1+c2)./2;
  %c(:,k+numel(refindx)) = c2;
end


function [M] = ft_connectivity_mim(inputdata, varargin)

% FT_CONNECTIVITY_MIM computes the multivariate interaction measure from a
% data-matrix containing the cross-spectral density. This implements the method
% described in Ewald et al., Estimating true brain connectivity from EEG/MEG data
% invariant to linear and static trasformations in sensor space. Neuroimage, 2012;
% 476:488.
%
% Use as
%   [m] = hcp_connectivity_mim(inputdata, ...)
%
% The input data should be an array organized as
%   Channel x Channel x Frequency
%
% The output m contains the newChannel x newChannel x Frequency connectivity measure,
% with newChannel equal to max(indices).
%
% Additional optional input arguments come as key-value pairs:
%   'indices' = 1xN vector with indices of the groups to which the channels belong,
%               e.g. [1 1 2 2] for a 2-by-2 connectivity between 2 planar MEG channels.
%
%
% See also CONNECTIVITY, FT_CONNECTIVITYANALYSIS

% Copyright (C) 2011-2014 by the Human Connectome Project, WU-Minn Consortium (1U54MH091657)
% Copyright (C) 2021 Jan-Mathijs Schoffelen, DCCN
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

indices = ft_getopt(varargin, 'indices');

if isempty(indices) && isequal(size(inputdata(:,:,1)), [2 2])
  % simply assume two channels
  indices = [1 1 2 2];
end

sizein  = size(inputdata);
sizeout = [sizein 1];
sizeout(1:2) = max(indices);

% compute the inverse of the auto terms only once for speed up
invC = cell(sizeout(1),sizeout(3));
for kk = 1:sizeout(3)
  for k = 1:sizeout(1)
    invC{k,kk} = pinv(real(inputdata(indices==k,indices==k,kk)));
  end
end

M = zeros(sizeout);
for kk = 1:sizeout(3)
  for k = 1:sizeout(1)
    for m = 1:sizeout(1)
      indx1 = indices==k;
      indx2 = indices==m;
      %cs_aa_re = real(input(indx1,indx1));
      %cs_bb_re = real(input(indx2,indx2));
      cs_ab_im = imag(inputdata(indx1,indx2,kk));
      
      %inv_cs_bb_re = pinv(cs_bb_re);
      %inv_cs_aa_re = pinv(cs_aa_re);
      inv_cs_bb_re = invC{m,kk};
      inv_cs_aa_re = invC{k,kk};
      transp_cs_ab_im = transpose(cs_ab_im);
      M(k,m,kk) = trace(inv_cs_aa_re*cs_ab_im*inv_cs_bb_re*transp_cs_ab_im); % try to speed up by dividing calculation in steps
    end
  end
end

% taking the mldivide and mrdivide operators doesn't change the results, but speeds up by a factor of 4 over 1000 iterations (on LM Notebook)
% m = trace(cs_aa_re\cs_ab_im*inv_cs_bb_re*transp_cs_ab_im);

% % % block_a=cs_aa_re\cs_ab_im;
% % % block_b=cs_bb_re\transpose(cs_ab_im);
% % % m = trace(block_a*block_b);

% m = trace(cs_aa_re\cs_ab_im*(cs_bb_re\transpose(cs_ab_im)));

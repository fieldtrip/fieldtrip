function [pow, csd, fourier] = timelock2freq(mom)

% TIMELOCK2FREQ transform the reconstructed dipole moment into
% something that again resembles the physical input parameter in
% the frequency domain. 
%
% This is needed after source reconstruction using FREQ2TIMELOCK.

% Copyright (C) 2005, Robert Oostenveld
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

n = size(mom,2)/2;

re = mom(:,1:n);
im = mom(:,(n+1):end);
% reconstruct the source level complex fourier representation
fourier = re + i*im;
% construct the source level csd matrix 
csd = fourier * ctranspose(fourier);
% estimate the source power
pow = trace(csd);


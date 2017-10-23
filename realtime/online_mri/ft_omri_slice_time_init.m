function STM = ft_omri_slice_time_init(X0, TR, deltaT)

% function STM = ft_omri_slice_time_init(X0, TR, deltaT);
%
% Initialize simple slice time correction structure.
% The algorithm will use plain linear interpolation.

% Copyright (C) 2010, Stefan Klanke
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

STM.dims = size(X0);
if nargin == 3
  if STM.dims(3)~=length(deltaT) || any(deltaT>TR) || any(deltaT<0) 
    error 'You must pass in a correct timing description vector';
  end
else
  deltaT = (0:(STM.dims(3)-1))*TR/STM.dims(3);
end

STM.deltaT = deltaT(:)';
STM.TR = TR;

xy = STM.dims(1)*STM.dims(2);

STM.weightOld = reshape(ones(xy,1)*(deltaT/TR), STM.dims);
STM.weightNew = reshape(ones(xy,1)*(1-deltaT/TR), STM.dims);

STM.lastScan = X0;

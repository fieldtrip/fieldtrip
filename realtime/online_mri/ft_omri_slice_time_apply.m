function [STM, Xs] = ft_omri_slice_time_apply(STM, X)

% function [STM, Xs] = ft_omri_slice_time_apply(STM, X)
%
% Put new scan X through slice time correction, by linear interpolation
% with last scan. The return value Xs is the signal sampled at deltaT = 0
% relative to the most recent scan.

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

sz = size(X);
Xs = X.*STM.weightNew + STM.lastScan.*STM.weightOld;
STM.lastScan = X;

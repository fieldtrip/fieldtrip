function ft_omri_quality_plot(motEst, curScan, varScan, maxVal, maxVar)

% function ft_omri_quality_plot(motEst, curScan, varScan, [maxVal, maxVar])
%
% motEst  should be Nx6 matrix of estimated motion parameters
% curScan should be [MxNxS] volume
% varScan should be [MxNxS] representation of the variation of the scans

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

if nargin<4
	maxVal = max(curScan(:));
end

if nargin<5
	maxVar = max(varScan(:));
end

trans = motEst(:,1:3);
rot   = motEst(:,4:6);

A1 = subplot(2,2,1);
plot(trans);
ylabel('Translation (mm)');
title('blue = X/readout dir.  |  green = Y/phase dir.  |  red = Z/slice axis');

A2 = subplot(2,2,2);
plot(rot);
ylabel('Rotation (deg)');
xlabel('scan number');

A3 = subplot(2,2,3);
imagesc(ft_omri_volume_to_mosaic(curScan), [0 maxVal]);
colormap(gray);
axis equal
axis off
title('Current scan');

A4 = subplot(2,2,4);
imagesc(ft_omri_volume_to_mosaic(varScan), [0 maxVar]);
colormap(gray);
axis equal
axis off
title('Variation');

set(A1,'Position',[0.07 0.75 0.9 0.2])
set(A2,'Position',[0.07 0.50 0.9 0.2])

set(A3,'Position',[0.05 0.001 0.45 0.45])
set(A4,'Position',[0.55 0.001 0.45 0.45])

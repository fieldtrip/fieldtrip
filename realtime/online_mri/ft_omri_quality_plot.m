function ft_omri_quality_plot(motEst, curScan, varScan, maxVal, maxVar)
%function ft_omri_quality_plot(motEst, curScan, varScan, [maxVal, maxVar])
%
%motEst should be Nx6 matrix of estimated motion parameters
%curScan a [MxNxS] volume
%varScan a [MxNxS] representation of the variation of the scans

% (C) 2010 S. Klanke

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

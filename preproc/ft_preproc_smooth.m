function [datsmooth] = ft_preproc_smooth(dat, n)

% FT_PREPROC_SMOOTH performs boxcar smoothing with specified length.
% Edge behavior is improved by implicit padding with the mean over
% half the boxcar length at the edges of the data segment.
%
% Use as
%   datsmooth = ft_preproc_smooth(dat, n)
%
% Where dat is an Nchan x Ntimepoints data matrix, and n the length
% of the boxcar smoothing kernel

% deal with padding
pad       = ceil(n/2);
edgeleft  = mean(dat(:,1:pad),2);
edgeright = mean(dat(:,(end+1-pad):end),2);
dat       = [edgeleft*ones(1,pad) dat edgeright*ones(1,pad)];

% create smoothing kernel
krn       = ones(1,n)/n;

% do the smoothing
datsmooth = convn(dat, krn, 'same');

% cut the eges
datsmooth = datsmooth(:,(pad+1):(end-pad));


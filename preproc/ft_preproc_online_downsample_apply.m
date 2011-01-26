function [DM, xd] = ft_preproc_online_downsample_apply(DM, x)
% function [DM, xd] = ft_preproc_online_downsample_apply(DM, x)
%
% Passes signal x (channels times samples) through the downsampler.
% Returns updated downsample model (numSkip!) and downsampled signal.
%
% 2010 S. Klanke

[dimX, numX] = size(x);

N  = size(x,2);
% to get number K of sample we can write out, subtract the skipped samples, 
% and then add maximum possible number of skip samples for next time (=DM.factor-1)
K  = floor((N - DM.numSkip + DM.factor-1)/DM.factor);
		
startIdx = 1+DM.numSkip;
endIdx   = 1+DM.numSkip + (K-1)*DM.factor;

xd = x(:,startIdx:DM.factor:endIdx);
DM.numSkip = DM.factor-1-(N-endIdx); % for next time
	
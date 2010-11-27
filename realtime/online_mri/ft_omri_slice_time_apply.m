function [STM, Xs] = ft_omri_slice_time_apply(STM, X)
% function [STM, Xs] = ft_omri_slice_time_apply(STM, X)
%
% Put new scan X through slice time correction, by linear interpolation
% with last scan. The return value Xs is the signal sampled at deltaT = 0
% relative to the most recent scan.

% 2010 S.Klanke

sz = size(X);
Xs = X.*STM.weightNew + STM.lastScan.*STM.weightOld;
STM.lastScan = X;
function [pow, csd, fourier] = timelock2freq(mom);

% TIMELOCK2FREQ transform the reconstructed dipole moment into
% something that again resembles the physical input parameter in
% the frequency domain. 
%
% This is needed after source reconstruction using FREQ2TIMELOCK.

% Copyright (C) 2005, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

n = size(mom,2)/2;

re = mom(:,1:n);
im = mom(:,(n+1):end);
% reconstruct the source level complex fourier representation
fourier = re + i*im;
% construct the source level csd matrix 
csd = fourier * ctranspose(fourier);
% estimate the source power
pow = trace(csd);


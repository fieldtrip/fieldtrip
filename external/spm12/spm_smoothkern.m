function krn = spm_smoothkern(fwhm,x,t)
% Generate a Gaussian smoothing kernel
% FORMAT krn = spm_smoothkern(fwhm,x,t)
% fwhm - full width at half maximum
% x    - position
% t    - either 0 (nearest neighbour) or 1 (linear).
%        [Default: 1]
%
% krn  - value of kernel at position x
%__________________________________________________________________________
%
% For smoothing images, one should really convolve a Gaussian with a sinc
% function. For smoothing histograms, the kernel should be a Gaussian
% convolved with the histogram basis function used. This function returns
% a Gaussian convolved with a triangular (1st degree B-spline) basis 
% function (by default). A Gaussian convolved with a hat function (0th 
% degree B-spline) can also be returned.
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_smoothkern.m 4419 2011-08-03 18:42:35Z guillaume $


if nargin<3, t = 1; end

% Variance from FWHM
s = (fwhm/sqrt(8*log(2)))^2+eps;

% The simple way to do it. Not good for small FWHM
% krn = (1/sqrt(2*pi*s))*exp(-(x.^2)/(2*s));

if t==0
    % Gaussian convolved with 0th degree B-spline
    % int(exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t= -0.5..0.5)
    w1  = 1/sqrt(2*s);
    krn = 0.5*(erf(w1*(x+0.5))-erf(w1*(x-0.5)));
    krn(krn<0) = 0;

elseif t==1
    % Gaussian convolved with 1st degree B-spline
    %  int((1-t)*exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t= 0..1)
    % +int((t+1)*exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t=-1..0)
    w1  =  0.5*sqrt(2/s);
    w2  = -0.5/s;
    w3  = sqrt(s/2/pi);
    krn = 0.5*(erf(w1*(x+1)).*(x+1) + erf(w1*(x-1)).*(x-1) - 2*erf(w1*x   ).* x)...
          +w3*(exp(w2*(x+1).^2)     + exp(w2*(x-1).^2)     - 2*exp(w2*x.^2));
    krn(krn<0) = 0;

else
    error('Only defined for nearest neighbour and linear interpolation.');
    % If anyone knows a nice formula for a sinc function convolved with a
    % a Gaussian, then that could be quite useful.
end

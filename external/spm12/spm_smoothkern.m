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
%
% Note that the convolution kernel returned by this function differ from
% the ones that other packages currently use for Gaussian smoothing -
% particularly when the FWHM is small compared with the voxel dimensions.
% The fact that SPM does it differently from other software does not mean
% that it is wrong.
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


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
    % This should probably be based on https://arxiv.org/pdf/1608.05854.pdf:
    %     Martin TB, Prunet S, Drissen L. Optimal fitting of Gaussian-apodized
    %     or under-resolved emission lines in Fourier transform spectra
    %     providing new insights on the velocity structure of NGC 6720. Monthly
    %     Notices of the Royal Astronomical Society. 2016 Sep 14;463(4):4223-38.
    % (thanks for the pointer Guillaume).
end

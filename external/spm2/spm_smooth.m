function spm_smooth(P,Q,s)
% 3 dimensional convolution of an image
% FORMAT spm_smooth(P,Q,S)
% P  - image to be smoothed
% Q  - filename for smoothed image
% S  - [sx sy sz] Guassian filter width {FWHM} in mm
%____________________________________________________________________________
%
% spm_smooth is used to smooth or convolve images in a file (maybe).
%
% The sum of kernel coeficients are set to unity.  Boundary
% conditions assume data does not exist outside the image in z (i.e.
% the kernel is truncated in z at the boundaries of the image space). S
% can be a vector of 3 FWHM values that specifiy an anisotropic
% smoothing.  If S is a scalar isotropic smoothing is implemented.
%
% If Q is not a string, it is used as the destination of the smoothed
% image.  It must already be defined with the same number of elements 
% as the image. 
%
%_______________________________________________________________________
% @(#)spm_smooth.m	2.3 John Ashburner, Tom Nichols 02/08/13

%-----------------------------------------------------------------------
if length(s) == 1; s = [s s s]; end

if ischar(P), P = spm_vol(P); end;
if isstruct(P),
	VOX = sqrt(sum(P.mat(1:3,1:3).^2));
else,
	VOX = [1 1 1];
end;

if ischar(Q) & isstruct(P),
	q         = Q;
	Q         = P;
	Q.fname   = q;
	if ~isfield(Q,'descrip'), Q.descrip = sprintf('SPM compatible'); end;
	Q.descrip = sprintf('%s - conv(%g,%g,%g)',Q.descrip, s);
end;

% compute parameters for spm_conv_vol
%-----------------------------------------------------------------------
s  = s./VOX;					% voxel anisotropy
s  = max(s,ones(size(s)));			% lower bound on FWHM
s  = s/sqrt(8*log(2));				% FWHM -> Gaussian parameter

x  = round(6*s(1)); x = [-x:x];
y  = round(6*s(2)); y = [-y:y];
z  = round(6*s(3)); z = [-z:z];
x  = exp(-(x).^2/(2*(s(1)).^2)); 
y  = exp(-(y).^2/(2*(s(2)).^2)); 
z  = exp(-(z).^2/(2*(s(3)).^2));
x  = x/sum(x);
y  = y/sum(y);
z  = z/sum(z);

i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;


if isstruct(Q), Q = spm_create_vol(Q); end; 
spm_conv_vol(P,Q,x,y,z,-[i,j,k]);
if isstruct(Q), Q = spm_close_vol(Q);  end;

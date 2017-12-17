function spm_smooth(P,Q,s,dtype)
% 3 dimensional convolution of an image
% FORMAT spm_smooth(P,Q,s,dtype)
% P     - image(s) to be smoothed (or 3D array)
% Q     - filename for smoothed image (or 3D array)
% s     - [sx sy sz] Gaussian filter width {FWHM} in mm (or edges)
% dtype - datatype [Default: 0 == same datatype as P]
%__________________________________________________________________________
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
%__________________________________________________________________________
% Copyright (C) 1994-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner & Tom Nichols
% $Id: spm_smooth.m 4419 2011-08-03 18:42:35Z guillaume $


%-Parameters & Arguments
%--------------------------------------------------------------------------
if numel(s) == 1, s = [s s s]; end
if nargin < 4, dtype = 0; end

if ischar(P), P = spm_vol(P); end

if isstruct(P)
    for i= 1:numel(P)
        smooth1(P(i),Q,s,dtype);
    end
else
    smooth1(P,Q,s,dtype);
end


%==========================================================================
% function smooth1(P,Q,s,dtype)
%==========================================================================
function smooth1(P,Q,s,dtype)

if isstruct(P)
    VOX = sqrt(sum(P.mat(1:3,1:3).^2));
else
    VOX = [1 1 1];
end

if ischar(Q) && isstruct(P)
    [pth,nam,ext,num] = spm_fileparts(Q);
    Q         = P;
    Q.fname   = fullfile(pth,[nam,ext]);
    if ~isempty(num)
        Q.n   = str2num(num);
    end
    if ~isfield(Q,'descrip'), Q.descrip = sprintf('SPM compatible'); end
    Q.descrip = sprintf('%s - conv(%g,%g,%g)',Q.descrip, s);

    if dtype~=0 % Need to figure out some rescaling
        Q.dt(1) = dtype;
        if ~isfinite(spm_type(Q.dt(1),'maxval')),
            Q.pinfo = [1 0 0]'; % float or double, so scalefactor of 1
        else
            % Need to determine the range of intensities
            if isfinite(spm_type(P.dt(1),'maxval')),
                % Integer types have a defined maximum value
                maxv = spm_type(P.dt(1),'maxval')*P.pinfo(1) + P.pinfo(2);
            else
                % Need to find the max and min values in original image
                mx = -Inf;
                mn =  Inf;
                for pl=1:P.dim(3)
                    tmp = spm_slice_vol(P,spm_matrix([0 0 pl]),P.dim(1:2),0);
                    tmp = tmp(isfinite(tmp));
                    mx  = max(max(tmp),mx);
                    mn  = min(min(tmp),mn);
                end
                maxv = max(mx,-mn);
            end
            sf      = maxv/spm_type(Q.dt(1),'maxval');
            Q.pinfo = [sf 0 0]';
        end
    end
end

%-Compute parameters for spm_conv_vol
%--------------------------------------------------------------------------
s  = s./VOX;                        % voxel anisotropy
s1 = s/sqrt(8*log(2));              % FWHM -> Gaussian parameter

x  = round(6*s1(1)); x = -x:x; x = spm_smoothkern(s(1),x,1); x  = x/sum(x);
y  = round(6*s1(2)); y = -y:y; y = spm_smoothkern(s(2),y,1); y  = y/sum(y);
z  = round(6*s1(3)); z = -z:z; z = spm_smoothkern(s(3),z,1); z  = z/sum(z);

i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;

if isstruct(Q), Q = spm_create_vol(Q); end

spm_conv_vol(P,Q,x,y,z,-[i,j,k]);

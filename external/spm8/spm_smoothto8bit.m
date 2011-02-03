function VO = spm_smoothto8bit(V,fwhm)
% 3 dimensional convolution of an image to 8bit data in memory
% FORMAT VO = spm_smoothto8bit(V,fwhm)
% V     - mapped image to be smoothed
% fwhm  - FWHM of Guassian filter width in mm
% VO    - smoothed volume in a form that can be used by the
%         spm_*_vol.mex* functions.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$


if nargin>1 & fwhm>0,
    VO = smoothto8bit(V,fwhm);
else,
    VO = V;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function VO = smoothto8bit(V,fwhm)
% 3 dimensional convolution of an image to 8bit data in memory
% FORMAT VO = smoothto8bit(V,fwhm)
% V     - mapped image to be smoothed
% fwhm  - FWHM of Guassian filter width in mm
% VO    - smoothed volume in a form that can be used by the
%         spm_*_vol.mex* functions.
%_______________________________________________________________________

vx   = sqrt(sum(V.mat(1:3,1:3).^2));
s    = (fwhm./vx./sqrt(8*log(2)) + eps).^2;
r    = cell(1,3);
for i=1:3,
    r{i}.s = ceil(3.5*sqrt(s(i)));
    x      = -r{i}.s:r{i}.s;
    r{i}.k = exp(-0.5 * (x.*x)/s(i))/sqrt(2*pi*s(i));
    r{i}.k = r{i}.k/sum(r{i}.k);
end;

buff = zeros([V.dim(1:2) r{3}.s*2+1]);

VO        = V;
VO.dt     = [spm_type('uint8') spm_platform('bigend')];
V0.dat    = uint8(0);
V0.dat(VO.dim(1:3)) = uint8(0);
VO.pinfo  = [];

for i=1:V.dim(3)+r{3}.s,
    if i<=V.dim(3),
        img      = spm_slice_vol(V,spm_matrix([0 0 i]),V.dim(1:2),0);
        msk      = find(~isfinite(img));
        img(msk) = 0;
        buff(:,:,rem(i-1,r{3}.s*2+1)+1) = ...
            conv2(conv2(img,r{1}.k,'same'),r{2}.k','same');
    else,
        buff(:,:,rem(i-1,r{3}.s*2+1)+1) = 0;
    end;

    if i>r{3}.s,
        kern    = zeros(size(r{3}.k'));
        kern(rem((i:(i+r{3}.s*2))',r{3}.s*2+1)+1) = r{3}.k';
        img     = reshape(buff,[prod(V.dim(1:2)) r{3}.s*2+1])*kern;
        img     = reshape(img,V.dim(1:2));
        ii      = i-r{3}.s;
        mx      = max(img(:));
        mn      = min(img(:));
        if mx==mn, mx=mn+eps; end;
        VO.pinfo(1:2,ii) = [(mx-mn)/255 mn]';
        VO.dat(:,:,ii)   = uint8(round((img-mn)*(255/(mx-mn))));
    end;
end;


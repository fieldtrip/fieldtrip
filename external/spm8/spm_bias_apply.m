function VO = spm_bias_apply(V,T)
% Apply bias field to image.
%
% FORMAT VO = spm_bias_apply(V,T)
%   V - filename or vol struct of image
%   T - DCT of bias field or filename containing this
%   VO - bias corrected volume structure.
%
% If no output arguments, then the bias corrected image is written to
% disk, prefixed by 'm'.
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_bias_apply.m 3756 2010-03-05 18:43:37Z guillaume $


if ischar(V),
    V = spm_vol(V);
end;

if ischar(T),
    s = load(T);
    T = s.T;
end;
nbas           = [size(T) 1];
nbas           = nbas(1:3);
B1             = spm_dctmtx(V(1).dim(1),nbas(1));
B2             = spm_dctmtx(V(1).dim(2),nbas(2));
B3             = spm_dctmtx(V(1).dim(3),nbas(3));

VO             = V;
VO.dt          = [spm_type('float32') spm_platform('bigend')];

if nargout==0,
    [pth,nm,xt,vr] = spm_fileparts(deblank(V.fname));
    VO.fname       = fullfile(pth,['m' nm xt vr]);
    %VO.fname      = ['m' nm xt vr];
    VO.pinfo       = [1 0 0]';
    VO             = spm_create_vol(VO);
else
    VO.fname       = 'bias_corrected.img';
    VO.pinfo       = [1 0]';
    VO.dat(1,1,1)  = single(0);
    VO.dat(VO.dim(1),VO.dim(2),VO.dim(3)) = 0;
end;

for p=1:V.dim(3),
    M   = spm_matrix([0 0 p]);
    img = spm_slice_vol(V, M, V.dim(1:2), 1);
    t   = reshape(T,  nbas(1)*nbas(2), nbas(3));
    t   = reshape(t*B3(p,:)', nbas(1), nbas(2));
    img = img.*exp(B1*t*B2');
    if nargout==0,
        VO  = spm_write_plane(VO,img,p);
    else
        VO.dat(:,:,p) = single(img);
    end;
end;
return;
%=======================================================================

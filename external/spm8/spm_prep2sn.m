function [po,pin] = spm_prep2sn(p)
% Convert the output from spm_preproc into an sn.mat
% FORMAT [po,pin] = spm_prep2sn(p)
% p   - the results of spm_preproc
% po  - the output in a form that can be used by
%       spm_write_sn.
% pin - the inverse transform in a form that can be
%       used by spm_write_sn.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$


if ischar(p), p = load(p); end;

VG          = p.tpm;
VF          = p.image;
[Y1,Y2,Y3]  = create_def(p.Twarp,VF,VG(1),p.Affine);
[Y1,Y2,Y3]  =     spm_invdef(Y1,Y2,Y3,VG(1).dim(1:3),eye(4),eye(4));
MT          =     procrustes(Y1,Y2,Y3,VG(1),VF.mat);
d           = size(p.Twarp);
[Affine,Tr] = reparameterise(Y1,Y2,Y3,VG(1),VF.mat,MT,max(d(1:3)+2,[8 8 8]));
flags       = struct(...
    'ngaus',    p.ngaus,...
    'mg',       p.mg,...
    'mn',       p.mn,...
    'vr',       p.vr,...
    'warpreg',  p.warpreg,...
    'warpco',   p.warpco,...
    'biasreg',  p.biasreg,...
    'biasfwhm', p.biasfwhm,...
    'regtype',  p.regtype,...
    'fudge',    p.fudge,...
    'samp',     p.samp,...
    'msk',      p.msk,...
    'Affine',   p.Affine,...
    'Twarp',    p.Twarp,...
    'Tbias',    p.Tbias,...
    'thresh',   p.thresh);

if nargout==0,
    [pth,nam,ext] = fileparts(VF.fname);
    fnam          = fullfile(pth,[nam '_seg_sn.mat']);
    if spm_matlab_version_chk('7') >= 0,
        save(fnam,'-V6','VG','VF','Tr','Affine','flags');
    else
        save(fnam,'VG','VF','Tr','Affine','flags');
    end;
else
    po = struct(...
        'VG',     VG,...
        'VF',     VF,...
        'Tr',     Tr,...
        'Affine', Affine,...
        'flags',  flags);
end;
if nargout>=2,
   % Parameterisation for the inverse
    pin = struct(...
        'VG',     p.image,...
        'VF',     p.tpm(1),...
        'Tr',     p.Twarp,...
        'Affine', p.tpm(1).mat\p.Affine*p.image.mat,...
        'flags',  flags);
    % VG = p.image;
    % VF = p.tpm(1);
    % Tr = p.Twarp;
    % Affine = p.tpm(1).mat\p.Affine*p.image.mat;
    % save('junk_sn.mat','VG','VF','Tr','Affine');
end;
return;
%=======================================================================

%=======================================================================
function [Affine,Tr] = reparameterise(Y1,Y2,Y3,B,M2,MT,d2)
% Take a deformation field and reparameterise in the same form
% as used by the spatial normalisation routines of SPM
d          = [size(Y1) 1];
[x1,x2,o]  = ndgrid(1:d(1),1:d(2),1);
x3         = 1:d(3);
Affine     = M2\MT*B(1).mat;
A          = inv(Affine);

B1  = spm_dctmtx(d(1),d2(1));
B2  = spm_dctmtx(d(2),d2(2));
B3  = spm_dctmtx(d(3),d2(3));
pd  = prod(d2(1:3));
AA  = eye(pd)*0.01;
Ab  = zeros(pd,3);
spm_progress_bar('init',length(x3),['Reparameterising'],'Planes completed');
mx = [0 0 0];
for z=1:length(x3),
    y1       = double(Y1(:,:,z));
    y2       = double(Y2(:,:,z));
    y3       = double(Y3(:,:,z));
    msk      = isfinite(y1);
    w        = double(msk);
    y1(~msk) = 0;
    y2(~msk) = 0;
    y3(~msk) = 0;
    z1       = A(1,1)*y1+A(1,2)*y2+A(1,3)*y3 + w.*(A(1,4) - x1);
    z2       = A(2,1)*y1+A(2,2)*y2+A(2,3)*y3 + w.*(A(2,4) - x2);
    z3       = A(3,1)*y1+A(3,2)*y2+A(3,3)*y3 + w *(A(3,4) - z );
    b3       = B3(z,:)';
    Ab(:,1)  = Ab(:,1) + kron(b3,spm_krutil(z1,B1,B2,0));
    Ab(:,2)  = Ab(:,2) + kron(b3,spm_krutil(z2,B1,B2,0));
    Ab(:,3)  = Ab(:,3) + kron(b3,spm_krutil(z3,B1,B2,0));
    AA       = AA  + kron(b3*b3',spm_krutil(w, B1,B2,1));
    spm_progress_bar('set',z);
end;
spm_progress_bar('clear');
Tr  = reshape(AA\Ab,[d2(1:3) 3]); drawnow;
return;
%=======================================================================

%=======================================================================
function MT = procrustes(Y1,Y2,Y3,B,M2)
% Take a deformation field and determine the closest rigid-body
% transform to match it, with weighing.
%
% Example Reference:
% F. L. Bookstein (1997).  "Landmark Methods for Forms Without
% Landmarks: Morphometrics of Group Differences in Outline Shape"
% Medical Image Analysis 1(3):225-243

M1         = B.mat;
d          = B.dim(1:3);
[x1,x2,o]  = ndgrid(1:d(1),1:d(2),1);
x3         = 1:d(3);
c1  = [0 0 0];
c2  = [0 0 0];
sw  =  0;
spm_progress_bar('init',length(x3),['Procrustes (1)'],'Planes completed');
for z=1:length(x3),
    y1         = double(Y1(:,:,z));
    y2         = double(Y2(:,:,z));
    y3         = double(Y3(:,:,z));
    msk        = find(isfinite(y1));
    w          = spm_sample_vol(B(1),x1(msk),x2(msk),o(msk)*z,0);
    swz        = sum(w(:));
    sw         = sw+swz;
    c1         = c1 + [w'*[x1(msk) x2(msk)] swz*z ];
    c2         = c2 +  w'*[y1(msk) y2(msk) y3(msk)];
    spm_progress_bar('set',z);
end;
spm_progress_bar('clear');
c1 = c1/sw;
c2 = c2/sw;
T1 = [eye(4,3) M1*[c1 1]'];
T2 = [eye(4,3) M2*[c2 1]'];
C  = zeros(3);
spm_progress_bar('init',length(x3),['Procrustes (2)'],'Planes completed');
for z=1:length(x3),
    y1         = double(Y1(:,:,z));
    y2         = double(Y2(:,:,z));
    y3         = double(Y3(:,:,z));
    msk        = find(isfinite(y1));
    w          = spm_sample_vol(B(1),x1(msk),x2(msk),o(msk)*z,0);
    C = C + [(x1(msk)-c1(1)).*w (x2(msk)-c1(2)).*w (    z-c1(3))*w ]' * ...
            [(y1(msk)-c2(1))    (y2(msk)-c2(2))    (y3(msk)-c2(3)) ];
    spm_progress_bar('set',z);
end;
spm_progress_bar('clear');
[u,s,v]    = svd(M1(1:3,1:3)*C*M2(1:3,1:3)');
R          = eye(4);
R(1:3,1:3) = v*u';
MT         = T2*R*inv(T1);
return;
%=======================================================================

%=======================================================================
function [Y1,Y2,Y3] = create_def(T,VG,VF,Affine)
% Generate a deformation field from its parameterisation.
d2   = size(T);
d    = VG.dim(1:3);
M    = VF.mat\Affine*VG.mat;
[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3   = 1:d(3);
B1   = spm_dctmtx(d(1),d2(1));
B2   = spm_dctmtx(d(2),d2(2));
B3   = spm_dctmtx(d(3),d2(3));
[pth,nam,ext] = fileparts(VG.fname);
spm_progress_bar('init',length(x3),['Creating Def: ' nam],'Planes completed');
for z=1:length(x3),
    [y1,y2,y3] = defs(T,z,B1,B2,B3,x1,x2,x3,M);
    Y1(:,:,z)  = single(y1);
    Y2(:,:,z)  = single(y2);
    Y3(:,:,z)  = single(y3);
    spm_progress_bar('set',z);
end;
spm_progress_bar('clear');
return;
%=======================================================================

%=======================================================================
function [x1,y1,z1] = defs(sol,z,B1,B2,B3,x0,y0,z0,M)
if ~isempty(sol),
    x1a = x0    + transf(B1,B2,B3(z,:),sol(:,:,:,1));
    y1a = y0    + transf(B1,B2,B3(z,:),sol(:,:,:,2));
    z1a = z0(z) + transf(B1,B2,B3(z,:),sol(:,:,:,3));
else
    x1a = x0;
    y1a = y0;
    z1a = z0;
end;
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
return;
%=======================================================================

%=======================================================================
function t = transf(B1,B2,B3,T)
d2 = [size(T) 1];
t1 = reshape(T, d2(1)*d2(2),d2(3)); drawnow;
t1 = reshape(t1*B3', d2(1), d2(2)); drawnow;
t  = B1*t1*B2';
return;
%=======================================================================

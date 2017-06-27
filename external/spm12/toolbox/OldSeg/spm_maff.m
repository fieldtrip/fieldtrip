function M = spm_maff(varargin)
% Affine registration to MNI space using mutual information
% FORMAT M = spm_maff(P,samp,x,b0,MF,M,regtyp,ff)
% P       - filename or structure handle of image
% x       - cell array of {x1,x2,x3}, where x1 and x2 are
%           co-ordinates (from ndgrid), and x3 is a list of
%           slice numbers to use
% b0      - a cell array of belonging probability images
%           (see spm_load_priors.m).
% MF      - voxel-to-world transform of belonging probability
%           images
% M       - starting estimates
% regtype - regularisation type
%           'mni'   - registration of European brains with MNI space
%           'eastern' - registration of East Asian brains with MNI space
%           'rigid' - rigid(ish)-body registration
%           'subj'  - inter-subject registration
%           'none'  - no regularisation
% ff      - a fudge factor (derived from the one above)
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_maff.m 4873 2012-08-30 19:06:26Z john $

[buf,MG] = loadbuf(varargin{1:2});
M        = affreg(buf, MG, varargin{2:end});

return;
%_______________________________________________________________________
%_______________________________________________________________________
function [buf,MG] = loadbuf(V,x)
if ischar(V), V = spm_vol(V); end;
x1 = x{1};
x2 = x{2};
x3 = x{3};
% Load the image
V         = spm_vol(V);
d         = V(1).dim(1:3);
o         = ones(size(x1));
d         = [size(x1) length(x3)];
g         = zeros(d);
spm_progress_bar('Init',V.dim(3),'Loading volume','Planes loaded');
for i=1:d(3)
    g(:,:,i) = spm_sample_vol(V,x1,x2,o*x3(i),0);
    spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');

% Convert the image to unsigned bytes
[mn,mx] = spm_minmax(g);
sw = warning('off','all');
for z=1:length(x3),
    gz         = g(:,:,z);
    buf(z).msk = gz>mn & isfinite(gz);
    buf(z).nm  = sum(buf(z).msk(:));
    gz         = double(gz(buf(z).msk));
    buf(z).g   = uint8(round(gz*(255/mx)));
end;
warning(sw);
MG = V.mat;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function [M,h0] = affreg(buf,MG,x,b0,MF,M,regtyp,ff)
% Do the work

x1 = x{1};
x2 = x{2};
x3 = x{3};
[mu,isig] = spm_affine_priors(regtyp);
mu   = [zeros(6,1) ; mu];
isig = [zeros(6,12) ; zeros(6,6) isig];
isig      =  isig*ff;
Alpha0    =  isig;
Beta0     = -isig*mu;

sol  = M2P(M);
sol1 = sol;
ll   = -Inf;
nsmp = sum(cat(1,buf.nm));
pr   = struct('b',[],'db1',[],'db2',[],'db3',[]);
spm_plot_convergence('Init','Registering','Log-likelihood','Iteration');

for iter=1:200
    penalty = (sol1-mu)'*isig*(sol1-mu);
    T       = MF\P2M(sol1)*MG;
    R       = derivs(MF,sol1,MG);
    y1a     = T(1,1)*x1 + T(1,2)*x2 + T(1,4);
    y2a     = T(2,1)*x1 + T(2,2)*x2 + T(2,4);
    y3a     = T(3,1)*x1 + T(3,2)*x2 + T(3,4);
    h0      = zeros(256,length(b0)-1)+eps;
    for i=1:length(x3),
        if ~buf(i).nm, continue; end;
        y1    = y1a(buf(i).msk) + T(1,3)*x3(i);
        y2    = y2a(buf(i).msk) + T(2,3)*x3(i);
        y3    = y3a(buf(i).msk) + T(3,3)*x3(i);
        for k=1:size(h0,2),
            pr(k).b = spm_sample_priors(b0{k},y1,y2,y3,k==length(b0));
            h0(:,k) = h0(:,k) + spm_hist(buf(i).g,pr(k).b);
        end;
    end;
    h1    = (h0+eps);
    ssh   = sum(h1(:));
    krn   = spm_smoothkern(2,(-256:256)',0);
    h1    = conv2(h1,krn,'same');
    h1    = h1/ssh;
    h2    = log2(h1./(sum(h1,2)*sum(h1,1)));
    ll1   = sum(sum(h0.*h2))/ssh - penalty/ssh;
    spm_plot_convergence('Set',ll1);
    if ll1-ll<1e-5, break; end;
    ll    = ll1;
    sol   = sol1;
    Alpha = zeros(12);
    Beta  = zeros(12,1);
    for i=1:length(x3),
        nz    = buf(i).nm;
        if ~nz, continue; end;
        msk   = buf(i).msk;
        gi    = double(buf(i).g)+1;
        y1    = y1a(msk) + T(1,3)*x3(i);
        y2    = y2a(msk) + T(2,3)*x3(i);
        y3    = y3a(msk) + T(3,3)*x3(i);

        dmi1  = zeros(nz,1);
        dmi2  = zeros(nz,1);
        dmi3  = zeros(nz,1);
        for k=1:size(h0,2),
            [pr(k).b, pr(k).db1, pr(k).db2, pr(k).db3] = spm_sample_priors(b0{k},y1,y2,y3,k==length(b0));
            tmp  = -h2(gi,k);
            dmi1 = dmi1 + tmp.*pr(k).db1;
            dmi2 = dmi2 + tmp.*pr(k).db2;
            dmi3 = dmi3 + tmp.*pr(k).db3;
        end;
        x1m = x1(msk);
        x2m = x2(msk);
        x3m = x3(i);
        A = [dmi1.*x1m dmi2.*x1m dmi3.*x1m...
             dmi1.*x2m dmi2.*x2m dmi3.*x2m...
             dmi1 *x3m dmi2 *x3m dmi3 *x3m...
             dmi1      dmi2      dmi3];
        Alpha = Alpha + A'*A;
        Beta  = Beta  + sum(A,1)';
    end;
    drawnow;
    Alpha = R'*Alpha*R;
    Beta  = R'*Beta;

    % Gauss-Newton update
    sol1  = (Alpha+Alpha0)\(Alpha*sol - Beta - Beta0);
end;

spm_plot_convergence('Clear');
M = P2M(sol);
return;
%_______________________________________________________________________
%_______________________________________________________________________
function P = M2P(M)
% Polar decomposition parameterisation of affine transform,
% based on matrix logs
J  = M(1:3,1:3);
V  = sqrtm(J*J');
R  = V\J;

lV = logm(V);
lR = -logm(R);
if sum(sum(imag(lR).^2))>1e-6
    error('Rotations by pi are still a problem.');
else
    lR = real(lR);
end
P       = zeros(12,1);
P(1:3)  = M(1:3,4);
P(4:6)  = lR([2 3 6]);
P(7:12) = lV([1 2 3 5 6 9]);
return;
%_______________________________________________________________________
%_______________________________________________________________________
function M = P2M(P)
% Polar decomposition parameterisation of affine transform,
% based on matrix logs

% Translations
D      = P(1:3);
D      = D(:);

% Rotation part
ind    = [2 3 6];
T      = zeros(3);
T(ind) = -P(4:6);
R      = expm(T-T');

% Symmetric part (zooms and shears)
ind    = [1 2 3 5 6 9];
T      = zeros(3);
T(ind) = P(7:12);
V      = expm(T+T'-diag(diag(T)));

M      = [V*R D ; 0 0 0 1];
return;
%_______________________________________________________________________
%_______________________________________________________________________
function R = derivs(MF,P,MG)
% Numerically compute derivatives of Affine transformation matrix w.r.t.
% changes in the parameters.
R = zeros(12,12);
M0 = MF\P2M(P)*MG;
M0 = M0(1:3,:);
for i=1:12
    dp     = 0.000000001;
    P1     = P;
    P1(i)  = P1(i) + dp;
    M1     = MF\P2M(P1)*MG;
    M1     = M1(1:3,:);
    R(:,i) = (M1(:)-M0(:))/dp;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function [h0,d1] = reg_unused(M)
% Try to analytically compute the first and second derivatives of a
% penalty function w.r.t. changes in parameters.  It works for first
% derivatives, but I couldn't make it work for the second derivs - so
% I gave up and tried a new strategy.

T   = M(1:3,1:3);
[U,S,V] = svd(T);
s   = diag(S);
h0  = sum(log(s).^2);
d1s = 2*log(s)./s;
%d2s = 2./s.^2-2*log(s)./s.^2;
d1  = zeros(12,1);

for j=1:3
    for i1=1:9
        T1     = zeros(3,3);
        T1(i1) = 1;
        t1     = U(:,j)'*T1*V(:,j);
        d1(i1) = d1(i1) + d1s(j)*t1;
    end;
end;
return;
%_______________________________________________________________________
%_______________________________________________________________________
function M = P2M_unused(P)
% SVD parameterisation of affine transform, based on matrix-logs.

% Translations
D      = P(1:3);
D      = D(:);

% Rotation U
ind    = [2 3 6];
T      = zeros(3);
T(ind) = P(4:6);
U      = expm(T-T');

% Diagonal zooming matrix
S      = expm(diag(P(7:9)));

% Rotation V'
T(ind) = P(10:12);
V      = expm(T'-T);

M      = [U*S*V' D ; 0 0 0 1];
return;
%_______________________________________________________________________
%_______________________________________________________________________


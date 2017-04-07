function [M,ll,h] = spm_maff8(varargin)
% Affine registration to MNI space using mutual information
% FORMAT [M,ll,h] = spm_maff8(P,samp,fwhm,tpm,M0,regtyp)
% P       - filename or structure handle of image
% samp    - distance between sample points (mm).  Small values are
%           better, but things run more slowly.
% fwhm    - smoothness estimate for computing a fudge factor.  Estimate
%           is a full width at half maximum of a Gaussian (in mm). 
% tpm     - data structure encoding a tissue probability map, generated
%           via spm_load_priors8.m.
% M0      - starting estimates for the affine transform (or [] to use
%           default values).
% regtype - regularisation type
%           'mni'     - registration of European brains with MNI space
%           'eastern' - registration of East Asian brains with MNI space
%           'rigid'   - rigid(ish)-body registration
%           'subj'    - inter-subject registration
%           'none'    - no regularisation
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_maff8.m 6421 2015-04-23 16:54:25Z john $

[buf,MG,x,ff] = loadbuf(varargin{1:3});
[M,ll,h]      = affreg(buf, MG, x, ff, varargin{4:end});

return;


%==========================================================================
% function [buf,MG,x,ff] = loadbuf(V,samp,fwhm)
%==========================================================================
function [buf,MG,x,ff] = loadbuf(V,samp,fwhm)
if ischar(V), V = spm_vol(V); end;
d       = V(1).dim(1:3);
vx      = sqrt(sum(V(1).mat(1:3,1:3).^2));  % Voxel sizes
sk      = max([1 1 1],round(samp*[1 1 1]./vx));
[x1,x2] = ndgrid(1:sk(1):d(1),1:sk(2):d(2));
x3      = 1:sk(3):d(3);

% Fudge Factor - to (approximately) account for
% non-independence of voxels
s    = (fwhm+mean(vx))/sqrt(8*log(2));                 % Standard deviation
ff   = prod(4*pi*(s./vx./sk).^2 + 1)^(1/2);

% Old version of fudge factor
%ff  = max(1,ff^3/prod(sk)/abs(det(V(1).mat(1:3,1:3))));

% Load the image
V         = spm_vol(V);
o         = ones(size(x1));
d         = [size(x1) length(x3)];
g         = zeros(d);
spm_progress_bar('Init',d(3),'Loading volume','Planes loaded');
for i=1:d(3)
    g(:,:,i) = spm_sample_vol(V,x1,x2,o*x3(i),0);
    spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');
spm_progress_bar('Init',d(3),'Initial Histogram','Planes complete');
mn = min(g(:));
mx = max(g(:));
sf = [mn 1;mx 1]\[1;4000];
h  = zeros(4000,1);
for i=1:d(3)
    p = g(:,:,i);
    p = p(isfinite(p) & (p~=0) & (p~=-3024));
    p = round(p*sf(1)+sf(2));
    h = h + accumarray(p(:),1,[4000,1]);
    spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');
spm_progress_bar('Init',d(3),'Converting to uint8','Planes complete');
h  = cumsum(h)/sum(h);
mn = (find(h>(0.0005),1)-sf(2))/sf(1);
mx = (find(h>(0.9995),1)-sf(2))/sf(1);
sf = [mn 1;mx 1]\[0;255];

if spm_type(V.dt(1),'intt'),
    scrand = V.pinfo(1);
    rand('seed',1);
else
    scrand = 0;
end

cl = cell(1,d(3));
buf = struct('nm',cl,'msk',cl,'g',cl);
for i=1:d(3),
    gz         = g(:,:,i);
    buf(i).msk = isfinite(gz) & gz~=0;
    buf(i).nm  = sum(buf(i).msk(:));
    if scrand, gz = gz + rand(size(gz))*scrand-scrand/2; end
    gz         = gz(buf(i).msk)*sf(1)+sf(2);
    buf(i).g   = uint8(max(min(round(gz),255),0));
    spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');
MG = V.mat;
x  = {x1,x2,x3};
return;


%==========================================================================
% function [M,h0] = affreg(buf,MG,x,ff,tpm,M,regtyp)
%==========================================================================
function [M,ll,h0] = affreg(buf,MG,x,ff,tpm,M,regtyp,maxit)
% Mutual Information Registration

if nargin<8, maxit=200; end

x1 = x{1};
x2 = x{2};
x3 = x{3};
[mu,isig] = spm_affine_priors(regtyp);
mu        = [zeros(6,1) ; mu];
Alpha0    = [eye(6,6)*0.00001 zeros(6,6) ; zeros(6,6) isig]*ff;

if ~isempty(M),
    sol  = M2P(M);
else
    sol  = mu;
end

ll   = -Inf;
krn  = spm_smoothkern(4,(-256:256)',0);

spm_plot_convergence('Init','Registering','Log-likelihood','Iteration');
h1       = ones(256,numel(tpm.dat));
nsearch  = 12;
for iter=1:200

    stepsize = 1;
    for search = 1:nsearch,
        if iter>1 
            sol1    = sol - stepsize*dsol;
        else
            sol1    = sol;
        end
        penalty = 0.5*(sol1-mu)'*Alpha0*(sol1-mu);
        M       = P2M(sol1);
        T       = tpm.M\M*MG;

       %fprintf('%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g | %g\n', sol1,penalty);
       %global st
       %st.vols{1}.premul = P2M(sol1);
       %spm_orthviews('Redraw')
       %drawnow

        y1a     = T(1,1)*x1 + T(1,2)*x2 + T(1,4);
        y2a     = T(2,1)*x1 + T(2,2)*x2 + T(2,4);
        y3a     = T(3,1)*x1 + T(3,2)*x2 + T(3,4);

        for i=1:length(x3)
            if ~buf(i).nm, continue; end
            y1    = y1a(buf(i).msk) + T(1,3)*x3(i);
            y2    = y2a(buf(i).msk) + T(2,3)*x3(i);
            y3    = y3a(buf(i).msk) + T(3,3)*x3(i);

            msk   = y3>=1;
            y1    = y1(msk);
            y2    = y2(msk);
            y3    = y3(msk);
            b     = spm_sample_priors8(tpm,y1,y2,y3);
            buf(i).b    = b;
            buf(i).msk1 = msk;
            buf(i).nm1  = sum(buf(i).msk1);
        end

        ll1 = 0;
        for subit=1:32
            h0  = zeros(256,numel(tpm.dat))+eps;
            if ~rem(subit,4),
                ll0 = ll1;
                ll1 = 0;
            end
            for i=1:length(x3),
                if ~buf(i).nm || ~buf(i).nm1, continue; end
                gm    = double(buf(i).g(buf(i).msk1))+1;
                q     = zeros(numel(gm),size(h0,2));
                for k=1:size(h0,2),
                    q(:,k) = h1(gm(:),k).*buf(i).b{k};
                end
                sq = sum(q,2) + eps;
                if ~rem(subit,4),
                    ll1 = ll1 + sum(log(sq));
                end
                for k=1:size(h0,2),
                    h0(:,k) = h0(:,k) + accumarray(gm,q(:,k)./sq,[256 1]);
                end
            end
            h1  = conv2((h0+eps)/sum(h0(:)),krn,'same');
            h1  = h1./(sum(h1,2)*sum(h1,1));
            if ~rem(subit,4),
                if (ll1-ll0)/sum(h0(:)) < 1e-5, break; end
            end
        end
        for i=1:length(x3)
            buf(i).b    = [];
            buf(i).msk1 = [];
        end

        ssh   = sum(h0(:));
        ll1   = (sum(sum(h0.*log(h1))) - penalty)/ssh/log(2);
       %fprintf('%g\t%g\n', sum(sum(h0.*log2(h1)))/ssh, -penalty/ssh);
       %spm_plot_convergence('Set',ll1);

        if iter==1, break; end               % No need for search
        if abs(ll1-ll)<1e-4, return; end     % Seems to have converged
        if (ll1<ll)
            stepsize = stepsize*0.5;        % Worse solution. Try again
            if search==nsearch, return; end; % Give up trying
        else
            break;                           % Better solution.  Carry on to next GN step.
        end
    end
    oh1   = h1;
    ll    = ll1;
    sol   = sol1;
    spm_plot_convergence('Set',ll);

    % Compute 1st and approximate second derivatives for computing Gauss-Newton update
    Alpha = zeros(12);   % 2nd derivatives with respect to an affine transform
    Beta  = zeros(12,1); % 1st derivatives with respect to an affine transform
    for i=1:length(x3)   % Loop over slices
        if ~buf(i).nm, continue; end
        gi    = double(buf(i).g)+1;
        y1    = y1a(buf(i).msk) + T(1,3)*x3(i);
        y2    = y2a(buf(i).msk) + T(2,3)*x3(i);
        y3    = y3a(buf(i).msk) + T(3,3)*x3(i);

        msk   = y3>=1;
        y1    = y1(msk);
        y2    = y2(msk);
        y3    = y3(msk);
        gi    = gi(msk);

        nz    = size(y1,1);
        if nz
            mi    = zeros(nz,1) + eps;
            dmi1  = zeros(nz,1);
            dmi2  = zeros(nz,1);
            dmi3  = zeros(nz,1);
            [b, db1, db2, db3] = spm_sample_priors8(tpm,y1,y2,y3);

            for k=1:size(h0,2)
                tmp  = h1(gi,k);
                mi   = mi   + tmp.*b{k};
                dmi1 = dmi1 + tmp.*db1{k};
                dmi2 = dmi2 + tmp.*db2{k};
                dmi3 = dmi3 + tmp.*db3{k};
            end
            dmi1 = dmi1./mi;
            dmi2 = dmi2./mi;
            dmi3 = dmi3./mi;

            % Convert from derivatives w.r.t. displacements at each voxel to
            % derivatives w.r.t. affine basis functions (at each voxel).
            x1m  = x1(buf(i).msk); x1m = x1m(msk);
            x2m  = x2(buf(i).msk); x2m = x2m(msk);
            x3m  = x3(i);
            A = [dmi1.*x1m dmi2.*x1m dmi3.*x1m...
                 dmi1.*x2m dmi2.*x2m dmi3.*x2m...
                 dmi1 *x3m dmi2 *x3m dmi3 *x3m...
                 dmi1      dmi2      dmi3];
            Alpha = Alpha + A'*A;      % Update Hessian
            Beta  = Beta  - sum(A,1)'; % Update gradient
        end
    end
    drawnow;

    % Convert from derivatives w.r.t. affine matrix, to derivatives w.r.t. parameters
    R     = derivs(tpm.M,sol,MG);
    Alpha = R'*Alpha*R;
    Beta  = R'*Beta;

    % Gauss-Newton increment direction
    dsol  = ((Alpha+Alpha0)\(Beta+Alpha0*(sol-mu)));
end

spm_plot_convergence('Clear');
M = P2M(sol);
return;


%==========================================================================
% function P = M2P(M)
%==========================================================================
function P = M2P(M)
% Polar decomposition parameterisation of affine transform,
% based on matrix logs
J  = M(1:3,1:3);
V  = sqrtm(J*J');
R  = V\J;

lV = logm(V);
lR = -logm(R);
if sum(sum(imag(lR).^2))>1e-6, error('Rotations by pi are still a problem.'); end;
P       = zeros(12,1);
P(1:3)  = M(1:3,4);
P(4:6)  = lR([2 3 6]);
P(7:12) = lV([1 2 3 5 6 9]);
P       = real(P);
return;


%==========================================================================
% function M = P2M(P)
%==========================================================================
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


%==========================================================================
% function R = derivs(MF,P,MG)
%==========================================================================
function R = derivs(MF,P,MG)
% Numerically compute derivatives of Affine transformation matrix w.r.t.
% changes in the parameters.
R  = zeros(12,12);
M0 = MF\P2M(P)*MG;
M0 = M0(1:3,:);
for i=1:12
    dp     = 0.0000001;
    P1     = P;
    P1(i)  = P1(i) + dp;
    M1     = MF\P2M(P1)*MG;
    M1     = M1(1:3,:);
    R(:,i) = (M1(:)-M0(:))/dp;
end
return;

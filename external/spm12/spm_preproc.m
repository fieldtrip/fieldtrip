function results = spm_preproc(varargin)
% Combined Segmentation and Spatial Normalisation
%
% FORMAT results = spm_preproc(V,opts)
%  V    - image to work with
%  opts - options
%  opts.tpm      - n tissue probability images for each class
%  opts.ngaus    - number of Gaussians per class (n+1 classes)
%  opts.warpreg  - warping regularisation
%  opts.warpco   - cutoff distance for DCT basis functions
%  opts.biasreg  - regularisation for bias correction
%  opts.biasfwhm - FWHM of Gausian form for bias regularisation
%  opts.regtype  - regularisation for affine part
%  opts.fudge    - a fudge factor
%  opts.msk      - unused
%__________________________________________________________________________
% Copyright (C) 2005-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_preproc.m 4916 2012-09-11 19:15:53Z guillaume $


SVNid     = '$Rev: 4916 $';

if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','OldSeg')); end
if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','OldNorm')); end

%-Say hello
%--------------------------------------------------------------------------
SPMid     = spm('FnBanner',mfilename,SVNid);

%-Parameters & Arguments
%--------------------------------------------------------------------------
opts0     = spm_get_defaults('old.preproc');
opts0     = rmfield(opts0,'output');
opts0.tpm = char(opts0.tpm); % In defaults, tpms are stored as cellstr
opts0.msk = '';

if ~nargin
    V = spm_select(1,'image');
else
    V = varargin{1};
end
if ischar(V), V = spm_vol(V); end

if nargin < 2
    opts = opts0;
else
    opts = varargin{2};
    fnms = fieldnames(opts0);
    for i=1:length(fnms)
        if ~isfield(opts,fnms{i}), opts.(fnms{i}) = opts0.(fnms{i}); end
    end
end

if length(opts.ngaus) ~= size(opts.tpm,1)+1
    error('Number of Gaussians per class is not compatible with number of classes');
end
K   = sum(opts.ngaus);
Kb  = length(opts.ngaus);
lkp = [];
for k=1:Kb
    lkp = [lkp ones(1,opts.ngaus(k))*k];
end

B         = spm_vol(opts.tpm);
b0        = spm_load_priors(B);

d         = V(1).dim(1:3);
vx        = sqrt(sum(V(1).mat(1:3,1:3).^2));
sk        = max([1 1 1],round(opts.samp*[1 1 1]./vx));
[x0,y0,o] = ndgrid(1:sk(1):d(1),1:sk(2):d(2),1);
z0        = 1:sk(3):d(3);
tiny      = eps;

vx   = sqrt(sum(V(1).mat(1:3,1:3).^2));
kron = inline('spm_krutil(a,b)','a','b');

% BENDING ENERGY REGULARIZATION for warping
%-----------------------------------------------------------------------
lam    = 0.001;
d2     = max(round((V(1).dim(1:3).*vx)/opts.warpco),[1 1 1]);
kx     = (pi*((1:d2(1))'-1)/d(1)/vx(1)).^2;
ky     = (pi*((1:d2(2))'-1)/d(2)/vx(2)).^2;
kz     = (pi*((1:d2(3))'-1)/d(3)/vx(3)).^2;
Cwarp  = (1*kron(kz.^2,kron(ky.^0,kx.^0)) +...
    1*kron(kz.^0,kron(ky.^2,kx.^0)) +...
    1*kron(kz.^0,kron(ky.^0,kx.^2)) +...
    2*kron(kz.^1,kron(ky.^1,kx.^0)) +...
    2*kron(kz.^1,kron(ky.^0,kx.^1)) +...
    2*kron(kz.^0,kron(ky.^1,kx.^1)) );
Cwarp  = Cwarp*opts.warpreg;
Cwarp  = [Cwarp*vx(1)^4 ; Cwarp*vx(2)^4 ; Cwarp*vx(3)^4];
Cwarp  = sparse(1:length(Cwarp),1:length(Cwarp),Cwarp,length(Cwarp),length(Cwarp));
B3warp = spm_dctmtx(d(3),d2(3),z0);
B2warp = spm_dctmtx(d(2),d2(2),y0(1,:)');
B1warp = spm_dctmtx(d(1),d2(1),x0(:,1));
lmR    = speye(size(Cwarp));
Twarp  = zeros([d2 3]);

% GAUSSIAN REGULARISATION for bias correction
%--------------------------------------------------------------------------
fwhm    = opts.biasfwhm;
sd      = vx(1)*V(1).dim(1)/fwhm; d3(1) = ceil(sd*2); krn_x   = exp(-(0:(d3(1)-1)).^2/sd.^2)/sqrt(vx(1));
sd      = vx(2)*V(1).dim(2)/fwhm; d3(2) = ceil(sd*2); krn_y   = exp(-(0:(d3(2)-1)).^2/sd.^2)/sqrt(vx(2));
sd      = vx(3)*V(1).dim(3)/fwhm; d3(3) = ceil(sd*2); krn_z   = exp(-(0:(d3(3)-1)).^2/sd.^2)/sqrt(vx(3));
Cbias   = kron(krn_z,kron(krn_y,krn_x)).^(-2)*opts.biasreg;
Cbias   = sparse(1:length(Cbias),1:length(Cbias),Cbias,length(Cbias),length(Cbias));
B3bias  = spm_dctmtx(d(3),d3(3),z0);
B2bias  = spm_dctmtx(d(2),d3(2),y0(1,:)');
B1bias  = spm_dctmtx(d(1),d3(1),x0(:,1));
lmRb    = speye(size(Cbias));
Tbias   = zeros(d3);


% Fudge Factor - to (approximately) account for non-independence of voxels
ff     = opts.fudge;
ff     = max(1,ff^3/prod(sk)/abs(det(V.mat(1:3,1:3))));
Cwarp  = Cwarp*ff;
Cbias  = Cbias*ff;

ll     = -Inf;
llr    = 0;
llrb   = 0;
tol1   = 1e-4; % Stopping criterion. For more accuracy, use a smaller value
d      = [size(x0) length(z0)];

f = zeros(d);
for z=1:length(z0)
    f(:,:,z) = spm_sample_vol(V,x0,y0,o*z0(z),0);
end
[thresh,mx]  = spm_minmax(f);
mn           = zeros(K,1);
% give same results each time
st = rand('state'); % st = rng;
rand('state',0); % rng(0,'v5uniform'); % rng('defaults');
for k1=1:Kb
    kk = sum(lkp==k1);
    mn(lkp==k1) = rand(kk,1)*mx;
end
rand('state',st); % rng(st);
vr           = ones(K,1)*mx^2;
mg           = ones(K,1)/K;

if ~isempty(opts.msk)
    VM = spm_vol(opts.msk);
    if sum(sum((VM.mat-V(1).mat).^2)) > 1e-6 || any(VM.dim(1:3) ~= V(1).dim(1:3))
        error('Mask must have the same dimensions and orientation as the image.');
    end
end

Affine  = eye(4);
if ~isempty(opts.regtype)
    Affine  = spm_maff(V,{x0,y0,z0},b0,B(1).mat,Affine,opts.regtype,ff*100);
    Affine  = spm_maff(V,{x0,y0,z0},b0,B(1).mat,Affine,opts.regtype,ff);
end
M       = B(1).mat\Affine*V(1).mat;

nm      = 0;
for z=1:length(z0)
    x1  = M(1,1)*x0 + M(1,2)*y0 + (M(1,3)*z0(z) + M(1,4));
    y1  = M(2,1)*x0 + M(2,2)*y0 + (M(2,3)*z0(z) + M(2,4));
    z1  = M(3,1)*x0 + M(3,2)*y0 + (M(3,3)*z0(z) + M(3,4));
    buf(z).msk = spm_sample_priors(b0{end},x1,y1,z1,1)<(1-1/512);
    fz         = f(:,:,z);
    %buf(z).msk = fz>thresh;
    buf(z).msk = buf(z).msk & isfinite(fz) & (fz~=0);
    
    if ~isempty(opts.msk),
        msk = spm_sample_vol(VM,x0,y0,o*z0(z),0);
        buf(z).msk = buf(z).msk & msk;
    end
    buf(z).nm  = sum(buf(z).msk(:));
    buf(z).f   = fz(buf(z).msk);
    nm         = nm + buf(z).nm;
    buf(z).bf(1:buf(z).nm,1) = single(1);
    buf(z).dat = single(0);
    if buf(z).nm,
        buf(z).dat(buf(z).nm,Kb) = single(0);
    end
end
clear f

finalit = 0;

spm_plot_convergence('Init','Processing','Log-likelihood','Iteration');
for iter=1:100
    
    if finalit
        % THIS CODE MAY BE USED IN FUTURE
        
        % Reload the data for the final iteration.  This iteration
        % does not do any registration, so there is no need to
        % mask out the background voxels.
        %------------------------------------------------------------------
        llrb  = -0.5*Tbias(:)'*Cbias*Tbias(:);
        for z=1:length(z0)
            fz         = spm_sample_vol(V,x0,y0,o*z0(z),0);
            buf(z).msk = fz~=0;
            if ~isempty(opts.msk)
                msk = spm_sample_vol(VM,x0,y0,o*z0(z),0);
                buf(z).msk = buf(z).msk & msk;
            end
            buf(z).nm  = sum(buf(z).msk(:));
            buf(z).f   = fz(buf(z).msk);
            nm         = nm + buf(z).nm;
            buf(z).bf(1:buf(z).nm,1) = single(1);
            buf(z).dat = single(0);
            if buf(z).nm
                buf(z).dat(buf(z).nm,Kb) = single(0);
            end
            if buf(z).nm
                bf        = transf(B1bias,B2bias,B3bias(z,:),Tbias);
                tmp       = bf(buf(z).msk);
                llrb      = llrb + sum(tmp);
                buf(z).bf = single(exp(tmp));
            end
        end
        
        % The background won't fit well any more, so increase the
        % variances of these Gaussians in order to give it a chance
        vr(lkp(K)) = vr(lkp(K))*8;
        
        spm_plot_convergence('Init','Processing','Log-likelihood','Iteration');
    end
    
    % Load the warped prior probability images into the buffer
    %----------------------------------------------------------------------
    for z=1:length(z0)
        if ~buf(z).nm, continue; end
        [x1,y1,z1] = defs(Twarp,z,B1warp,B2warp,B3warp,x0,y0,z0,M,buf(z).msk);
        for k1=1:Kb
            tmp              = spm_sample_priors(b0{k1},x1,y1,z1,k1==Kb);
            buf(z).dat(:,k1) = single(tmp);
        end
    end
    
    for iter1=1:10
        
        % Estimate cluster parameters
        %==================================================================
        for subit=1:40
            oll  = ll;
            mom0 = zeros(K,1)+tiny;
            mom1 = zeros(K,1);
            mom2 = zeros(K,1);
            mgm  = zeros(Kb,1);
            ll   = llr+llrb;
            for z=1:length(z0)
                if ~buf(z).nm, continue; end
                bf  = double(buf(z).bf);
                cr  = double(buf(z).f).*bf;
                q   =  zeros(buf(z).nm,K);
                b   =  zeros(buf(z).nm,Kb);
                s   =  zeros(buf(z).nm,1)+tiny;
                for k1=1:Kb
                    pr      = double(buf(z).dat(:,k1));
                    b(:,k1) = pr;
                    s       = s + pr*sum(mg(lkp==k1));
                end
                for k1=1:Kb
                    b(:,k1) = b(:,k1)./s;
                end
                mgm = mgm + sum(b,1)';
                for k=1:K
                    q(:,k) = mg(k)*b(:,lkp(k)) .* exp((cr-mn(k)).^2/(-2*vr(k)))/sqrt(2*pi*vr(k));
                end
                sq = sum(q,2)+tiny;
                ll = ll + sum(log(sq));
                for k=1:K % Moments
                    p1      = q(:,k)./sq; mom0(k) = mom0(k) + sum(p1(:));
                    p1      = p1.*cr;     mom1(k) = mom1(k) + sum(p1(:));
                    p1      = p1.*cr;     mom2(k) = mom2(k) + sum(p1(:));
                end
            end
            
            % Mixing proportions, Means and Variances
            for k=1:K
                mg(k) = (mom0(k)+eps)/(mgm(lkp(k))+eps);
                mn(k) = mom1(k)/(mom0(k)+eps);
                vr(k) =(mom2(k)-mom1(k)*mom1(k)/mom0(k)+1e6*eps)/(mom0(k)+eps);
                vr(k) = max(vr(k),eps);
            end
            
            if subit>1 || (iter>1 && ~finalit),
                spm_plot_convergence('Set',ll);
            end;
            if finalit, fprintf('Mix: %g\n',ll); end;
            if subit == 1
                ooll = ll;
            elseif (ll-oll)<tol1*nm
                % Improvement is small, so go to next step
                break;
            end
        end
        
        % Estimate bias
        %==================================================================
        if prod(d3)>0
            for subit=1:40
                % Compute objective function and its 1st and second derivatives
                Alpha = zeros(prod(d3),prod(d3)); % Second derivatives
                Beta  = zeros(prod(d3),1); % First derivatives
                ollrb = llrb;
                oll   = ll;
                ll    = llr+llrb;
                for z=1:length(z0)
                    if ~buf(z).nm, continue; end
                    bf  = double(buf(z).bf);
                    cr  = double(buf(z).f).*bf;
                    q   = zeros(buf(z).nm,K);
                    for k=1:K
                        q(:,k) = double(buf(z).dat(:,lkp(k)))*mg(k);
                    end
                    s = sum(q,2)+tiny;
                    for k=1:K
                        q(:,k) = q(:,k)./s .* exp((cr-mn(k)).^2/(-2*vr(k)))/sqrt(2*pi*vr(k));
                    end
                    sq    = sum(q,2)+tiny;
                    ll    = ll   + sum(log(sq));
                    
                    w1 = zeros(buf(z).nm,1);
                    w2 = zeros(buf(z).nm,1);
                    for k=1:K
                        tmp = q(:,k)./sq/vr(k);
                        w1  = w1 + tmp.*(mn(k) - cr);
                        w2  = w2 + tmp;
                    end
                    wt1   = zeros(d(1:2)); wt1(buf(z).msk) = 1 + cr.*w1;
                    wt2   = zeros(d(1:2)); wt2(buf(z).msk) = cr.*(cr.*w2 - w1);
                    b3    = B3bias(z,:)';
                    Beta  = Beta  + kron(b3,spm_krutil(wt1,B1bias,B2bias,0));
                    Alpha = Alpha + kron(b3*b3',spm_krutil(wt2,B1bias,B2bias,1));
                    clear w1 w2 wt1 wt2 b3
                end
                if finalit, fprintf('Bia: %g\n',ll); end
                if subit > 1 && ~(ll>oll)
                    % Hasn't improved, so go back to previous solution
                    Tbias = oTbias;
                    llrb  = ollrb;
                    for z=1:length(z0)
                        if ~buf(z).nm, continue; end
                        bf        = transf(B1bias,B2bias,B3bias(z,:),Tbias);
                        buf(z).bf = single(exp(bf(buf(z).msk)));
                    end
                    break;
                else
                    % Accept new solution
                    spm_plot_convergence('Set',ll);
                    oTbias = Tbias;
                    if subit > 1 && ~((ll-oll)>tol1*nm)
                        % Improvement is only small, so go to next step
                        break;
                    else
                        % Use new solution and continue the Levenberg-Marquardt iterations
                        Tbias  = reshape((Alpha + Cbias + lmRb)\((Alpha+lmRb)*Tbias(:) + Beta),d3);
                        llrb  = -0.5*Tbias(:)'*Cbias*Tbias(:);
                        for z=1:length(z0)
                            if ~buf(z).nm, continue; end
                            bf        = transf(B1bias,B2bias,B3bias(z,:),Tbias);
                            tmp       = bf(buf(z).msk);
                            llrb      = llrb + sum(tmp);
                            buf(z).bf = single(exp(tmp));
                        end
                    end
                end
            end
            if ~((ll-ooll)>tol1*nm), break; end
        end
    end
    
    if finalit, break; end
    
    % Estimate deformations
    %======================================================================
    mg1 = full(sparse(lkp,1,mg));
    ll  = llr+llrb;
    for z=1:length(z0)
        if ~buf(z).nm, continue; end
        bf  = double(buf(z).bf);
        cr  = double(buf(z).f).*bf;
        q   =  zeros(buf(z).nm,Kb);
        tmp =  zeros(buf(z).nm,1)+tiny;
        s   =  zeros(buf(z).nm,1)+tiny;
        for k1=1:Kb
            s = s + mg1(k1)*double(buf(z).dat(:,k1));
        end
        for k1=1:Kb
            kk = find(lkp==k1);
            pp = zeros(buf(z).nm,1);
            for k=kk
                pp = pp + exp((cr-mn(k)).^2/(-2*vr(k)))/sqrt(2*pi*vr(k))*mg(k);
            end
            q(:,k1) = pp;
            tmp     = tmp+pp.*double(buf(z).dat(:,k1))./s;
        end
        ll = ll + sum(log(tmp));
        for k1=1:Kb
            buf(z).dat(:,k1) = single(q(:,k1));
        end
    end
    
    for subit=1:20
        oll    = ll;
        A      = cell(3,3);
        A{1,1} = zeros(prod(d2));
        A{1,2} = zeros(prod(d2));
        A{1,3} = zeros(prod(d2));
        A{2,2} = zeros(prod(d2));
        A{2,3} = zeros(prod(d2));
        A{3,3} = zeros(prod(d2));
        Beta   = zeros(prod(d2)*3,1);
        
        for z=1:length(z0)
            if ~buf(z).nm, continue; end
            [x1,y1,z1] = defs(Twarp,z,B1warp,B2warp,B3warp,x0,y0,z0,M,buf(z).msk);
            b   = zeros(buf(z).nm,Kb);
            db1 = zeros(buf(z).nm,Kb);
            db2 = zeros(buf(z).nm,Kb);
            db3 = zeros(buf(z).nm,Kb);
            s   = zeros(buf(z).nm,1)+tiny;
            ds1 = zeros(buf(z).nm,1);
            ds2 = zeros(buf(z).nm,1);
            ds3 = zeros(buf(z).nm,1);
            p   = zeros(buf(z).nm,1)+tiny;
            dp1 = zeros(buf(z).nm,1);
            dp2 = zeros(buf(z).nm,1);
            dp3 = zeros(buf(z).nm,1);
            
            for k1=1:Kb
                [b(:,k1),db1(:,k1),db2(:,k1),db3(:,k1)] = spm_sample_priors(b0{k1},x1,y1,z1,k1==Kb);
                s   = s   + mg1(k1)*  b(:,k1);
                ds1 = ds1 + mg1(k1)*db1(:,k1);
                ds2 = ds2 + mg1(k1)*db2(:,k1);
                ds3 = ds3 + mg1(k1)*db3(:,k1);
            end
            for k1=1:Kb
                b(:,k1)   = b(:,k1)./s;
                db1(:,k1) = (db1(:,k1)-b(:,k1).*ds1)./s;
                db2(:,k1) = (db2(:,k1)-b(:,k1).*ds2)./s;
                db3(:,k1) = (db3(:,k1)-b(:,k1).*ds3)./s;
                
                pp  = double(buf(z).dat(:,k1));
                p   = p   + pp.*b(:,k1);
                dp1 = dp1 + pp.*(M(1,1)*db1(:,k1) + M(2,1)*db2(:,k1) + M(3,1)*db3(:,k1));
                dp2 = dp2 + pp.*(M(1,2)*db1(:,k1) + M(2,2)*db2(:,k1) + M(3,2)*db3(:,k1));
                dp3 = dp3 + pp.*(M(1,3)*db1(:,k1) + M(2,3)*db2(:,k1) + M(3,3)*db3(:,k1));
            end
            
            clear x1 y1 z1 b db1 db2 db3 s ds1 ds2 ds3
            
            tmp             = zeros(d(1:2));
            tmp(buf(z).msk) = dp1./p; dp1 = tmp;
            tmp(buf(z).msk) = dp2./p; dp2 = tmp;
            tmp(buf(z).msk) = dp3./p; dp3 = tmp;
            
            b3     = B3warp(z,:)';
            Beta   = Beta - [...
                kron(b3,spm_krutil(dp1,B1warp,B2warp,0))
                kron(b3,spm_krutil(dp2,B1warp,B2warp,0))
                kron(b3,spm_krutil(dp3,B1warp,B2warp,0))];
            
            b3b3   = b3*b3';
            A{1,1} = A{1,1} +  kron(b3b3,spm_krutil(dp1.*dp1,B1warp,B2warp,1));
            A{1,2} = A{1,2} +  kron(b3b3,spm_krutil(dp1.*dp2,B1warp,B2warp,1));
            A{1,3} = A{1,3} +  kron(b3b3,spm_krutil(dp1.*dp3,B1warp,B2warp,1));
            A{2,2} = A{2,2} +  kron(b3b3,spm_krutil(dp2.*dp2,B1warp,B2warp,1));
            A{2,3} = A{2,3} +  kron(b3b3,spm_krutil(dp2.*dp3,B1warp,B2warp,1));
            A{3,3} = A{3,3} +  kron(b3b3,spm_krutil(dp3.*dp3,B1warp,B2warp,1));
            
            clear b3 b3b3 tmp p dp1 dp2 dp3
        end
        
        Alpha = [A{1,1} A{1,2} A{1,3} ; A{1,2} A{2,2} A{2,3}; A{1,3} A{2,3} A{3,3}];
        clear A
        
        for subit1 = 1:3
            if iter==1,
                nTwarp = (Alpha+lmR*lam + 10*Cwarp)\((Alpha+lmR*lam)*Twarp(:) - Beta);
            else
                nTwarp = (Alpha+lmR*lam +    Cwarp)\((Alpha+lmR*lam)*Twarp(:) - Beta);
            end
            nTwarp = reshape(nTwarp,[d2 3]);
            nllr   = -0.5*nTwarp(:)'*Cwarp*nTwarp(:);
            nll    = nllr+llrb;
            for z=1:length(z0)
                if ~buf(z).nm, continue; end
                [x1,y1,z1] = defs(nTwarp,z,B1warp,B2warp,B3warp,x0,y0,z0,M,buf(z).msk);
                sq = zeros(buf(z).nm,1) + tiny;
                b  = zeros(buf(z).nm,Kb);
                s  = zeros(buf(z).nm,1)+tiny;
                for k1=1:Kb
                    b(:,k1) = spm_sample_priors(b0{k1},x1,y1,z1,k1==Kb);
                    s       = s + mg1(k1)*b(:,k1);
                end
                for k1=1:Kb
                    sq   = sq + double(buf(z).dat(:,k1)).*b(:,k1)./s;
                end
                clear b
                nll = nll + sum(log(sq));
                clear sq x1 y1 z1
            end
            if nll<ll
                % Worse solution, so use old solution and increase regularisation
                lam = lam*10;
            else
                % Accept new solution
                ll     = nll;
                llr    = nllr;
                Twarp  = nTwarp;
                lam    = lam*0.5;
                break
            end
        end
        
        spm_plot_convergence('Set',ll);
        if (ll-oll)<tol1*nm, break; end
    end
    
    if ~((ll-ooll)>tol1*nm)
        finalit = 1;
        break; % This can be commented out.
    end
end
spm_plot_convergence('Clear');

results        = opts;
results.image  = V;
results.tpm    = B;
results.Affine = Affine;
results.Twarp  = Twarp;
results.Tbias  = Tbias;
results.mg     = mg;
results.mn     = mn;
results.vr     = vr;
results.thresh = 0; %thresh;
results.ll     = ll;

fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#


%=======================================================================

%=======================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T),
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t  = zeros(size(B1,1),size(B2,1));
end;
return;
%=======================================================================

%=======================================================================
function [x1,y1,z1] = defs(Twarp,z,B1,B2,B3,x0,y0,z0,M,msk)
x1a = x0    + transf(B1,B2,B3(z,:),Twarp(:,:,:,1));
y1a = y0    + transf(B1,B2,B3(z,:),Twarp(:,:,:,2));
z1a = z0(z) + transf(B1,B2,B3(z,:),Twarp(:,:,:,3));
if nargin>=10,
    x1a = x1a(msk);
    y1a = y1a(msk);
    z1a = z1a(msk);
end;
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
return;
%=======================================================================

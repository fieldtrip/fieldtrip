function spm_preproc_write(p,opts)
% Write out VBM preprocessed data
% FORMAT spm_preproc_write(p,opts)
% p    - results from spm_prep2sn
% opts - writing options.  A struct containing these fields:
%        biascor - write bias corrected image
%        GM      - flags for which images should be written
%        WM      - similar to GM
%        CSF     - similar to GM
%____________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id$


if nargin==1,
    opts = struct('biascor',0,'GM',[0 0 1],'WM',[0 0 1],'CSF',[0 0 0],'cleanup',0);
end;
if numel(p)>0,
    b0  = spm_load_priors(p(1).VG);
end;
for i=1:numel(p),
    preproc_apply(p(i),opts,b0);
end;
return;
%=======================================================================

%=======================================================================
function preproc_apply(p,opts,b0)

%sopts = [opts.GM ; opts.WM ; opts.CSF];
nclasses = size(fieldnames(opts),1) - 2 ;
switch nclasses
    case 3
        sopts = [opts.GM ; opts.WM ; opts.CSF];
    case 4
       sopts = [opts.GM ; opts.WM ; opts.CSF ; opts.EXTRA1];
    case 5
       sopts = [opts.GM ; opts.WM ; opts.CSF ; opts.EXTRA1 ; opts.EXTRA2];
    otherwise
        error('######## unsupported number of classes....!!!')
end

[pth,nam,ext]=fileparts(p.VF.fname);
T    = p.flags.Twarp;
bsol = p.flags.Tbias;
d2   = [size(T) 1];
d    = p.VF.dim(1:3);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3  = 1:d(3);
d3  = [size(bsol) 1];
B1  = spm_dctmtx(d(1),d2(1));
B2  = spm_dctmtx(d(2),d2(2));
B3  = spm_dctmtx(d(3),d2(3));
bB3 = spm_dctmtx(d(3),d3(3),x3);
bB2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
bB1 = spm_dctmtx(d(1),d3(1),x1(:,1));

mg  = p.flags.mg;
mn  = p.flags.mn;
vr  = p.flags.vr;
K   = length(p.flags.mg);
Kb  = length(p.flags.ngaus);

for k1=1:size(sopts,1),
    %dat{k1}                = zeros(d(1:3),'uint8');
    dat{k1}                 = uint8(0);
    dat{k1}(d(1),d(2),d(3)) = 0;
    if sopts(k1,3),
        Vt        = struct('fname',   fullfile(pth,['c', num2str(k1), nam, ext]),...
                           'dim',     p.VF.dim,...
                           'dt',      [spm_type('uint8') spm_platform('bigend')],...
                           'pinfo',   [1/255 0 0]',...
                           'mat',     p.VF.mat,...
                           'n',       [1 1],...
                           'descrip', ['Tissue class ' num2str(k1)]);
        Vt        = spm_create_vol(Vt);
        VO(k1)    = Vt;
    end;
end;
if opts.biascor,
    VB = struct('fname',   fullfile(pth,['m', nam, ext]),...
                'dim',     p.VF.dim(1:3),...
                'dt',      [spm_type('float32') spm_platform('bigend')],...
                'pinfo',   [1 0 0]',...
                'mat',     p.VF.mat,...
                'n',       [1 1],...
                'descrip', 'Bias Corrected');
    VB = spm_create_vol(VB);
end;

lkp = []; for k=1:Kb, lkp = [lkp ones(1,p.flags.ngaus(k))*k]; end;

spm_progress_bar('init',length(x3),['Working on ' nam],'Planes completed');
M = p.VG(1).mat\p.flags.Affine*p.VF.mat;

for z=1:length(x3),

    % Bias corrected image
    f          = spm_sample_vol(p.VF,x1,x2,o*x3(z),0);
    cr         = exp(transf(bB1,bB2,bB3(z,:),bsol)).*f;
    if opts.biascor,
        % Write a plane of bias corrected data
        VB = spm_write_plane(VB,cr,z);
    end;

    if any(sopts(:)),
        msk        = (f==0) | ~isfinite(f);
        [t1,t2,t3] = defs(T,z,B1,B2,B3,x1,x2,x3,M);
        q          = zeros([d(1:2) Kb]);
        bt         = zeros([d(1:2) Kb]);
        for k1=1:Kb,
            bt(:,:,k1) = spm_sample_priors(b0{k1},t1,t2,t3,k1==Kb);
        end;
        b = zeros([d(1:2) K]);
        for k=1:K,
            b(:,:,k) = bt(:,:,lkp(k))*mg(k);
        end;
        s = sum(b,3);
        for k=1:K,
            p1            = exp((cr-mn(k)).^2/(-2*vr(k)))/sqrt(2*pi*vr(k)+eps);
            q(:,:,lkp(k)) = q(:,:,lkp(k)) + p1.*b(:,:,k)./s;
        end;
        sq = sum(q,3)+eps;
        sw = warning('off','MATLAB:divideByZero');
        for k1=1:size(sopts,1),
            tmp            = q(:,:,k1);
            tmp(msk)       = 0;
            tmp            = tmp./sq;
            dat{k1}(:,:,z) = uint8(round(255 * tmp));
        end;
        warning(sw);
    end;
    spm_progress_bar('set',z);
end;
spm_progress_bar('clear');

if opts.cleanup > 0,
    [dat{1},dat{2},dat{3}] = clean_gwc(dat{1},dat{2},dat{3}, opts.cleanup);
end;
if any(sopts(:,3)),
    for z=1:length(x3),
        for k1=1:size(sopts,1),
            if sopts(k1,3),
                tmp = double(dat{k1}(:,:,z))/255;
                spm_write_plane(VO(k1),tmp,z);
            end;
        end;
    end;
end;

for k1=1:size(sopts,1),
    if any(sopts(k1,1:2)),
        so      = struct('wrap',[0 0 0],'interp',1,'vox',[NaN NaN NaN],...
                         'bb',ones(2,3)*NaN,'preserve',0);
        ovx     = abs(det(p.VG(1).mat(1:3,1:3)))^(1/3);
        fwhm    = max(ovx./sqrt(sum(p.VF.mat(1:3,1:3).^2))-1,0.1);
        dat{k1} = decimate(dat{k1},fwhm);
        fn      = fullfile(pth,['c', num2str(k1), nam, ext]);
        dim     = [size(dat{k1}) 1];
        VT      = struct('fname',fn,'dim',dim(1:3),...
                     'dt', [spm_type('uint8') spm_platform('bigend')],...
                     'pinfo',[1/255 0]','mat',p.VF.mat,'dat',dat{k1});
        if sopts(k1,2),
            spm_write_sn(VT,p,so);
        end;
        so.preserve = 1;
        if sopts(k1,1),
           VN       = spm_write_sn(VT,p,so);
           VN.fname = fullfile(pth,['mwc', num2str(k1), nam, ext]);
           spm_write_vol(VN,VN.dat);
        end;
    end;
end;
return;
%=======================================================================

%=======================================================================
function [x1,y1,z1] = defs(sol,z,B1,B2,B3,x0,y0,z0,M)
x1a = x0    + transf(B1,B2,B3(z,:),sol(:,:,:,1));
y1a = y0    + transf(B1,B2,B3(z,:),sol(:,:,:,2));
z1a = z0(z) + transf(B1,B2,B3(z,:),sol(:,:,:,3));
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
return;
%=======================================================================

%=======================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t = zeros(size(B1,1),size(B2,1),size(B3,1));
end;
return;
%=======================================================================

%=======================================================================
function dat = decimate(dat,fwhm)
% Convolve the volume in memory (fwhm in voxels).
lim = ceil(2*fwhm);
x  = -lim(1):lim(1); x = spm_smoothkern(fwhm(1),x); x  = x/sum(x);
y  = -lim(2):lim(2); y = spm_smoothkern(fwhm(2),y); y  = y/sum(y);
z  = -lim(3):lim(3); z = spm_smoothkern(fwhm(3),z); z  = z/sum(z);
i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
spm_conv_vol(dat,dat,x,y,z,-[i j k]);
return;
%=======================================================================

%=======================================================================
function [g,w,c] = clean_gwc(g,w,c, level)
if nargin<4, level = 1; end;

b    = w;
b(1) = w(1);

% Build a 3x3x3 seperable smoothing kernel
%-----------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

th1 = 0.15;
if level==2, th1 = 0.2; end;
% Erosions and conditional dilations
%-----------------------------------------------------------------------
niter = 32;
spm_progress_bar('Init',niter,'Extracting Brain','Iterations completed');
for j=1:niter,
        if j>2, th=th1; else th=0.6; end; % Dilate after two its of erosion.
        for i=1:size(b,3),
                gp = double(g(:,:,i));
                wp = double(w(:,:,i));
                bp = double(b(:,:,i))/255;
                bp = (bp>th).*(wp+gp);
                b(:,:,i) = uint8(round(bp));
        end;
        spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
        spm_progress_bar('Set',j);
end;
th = 0.05;
for i=1:size(b,3),
        gp       = double(g(:,:,i))/255;
        wp       = double(w(:,:,i))/255;
        cp       = double(c(:,:,i))/255;
        bp       = double(b(:,:,i))/255;
        bp       = ((bp>th).*(wp+gp))>th;
        g(:,:,i) = uint8(round(255*gp.*bp./(gp+wp+cp+eps)));
        w(:,:,i) = uint8(round(255*wp.*bp./(gp+wp+cp+eps)));
        c(:,:,i) = uint8(round(255*(cp.*bp./(gp+wp+cp+eps)+cp.*(1-bp))));
end;
spm_progress_bar('Clear');
return;
%=======================================================================

%=======================================================================


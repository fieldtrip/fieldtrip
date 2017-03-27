function out = spm_dartel_import(job)
% Import subjects' data for use with Dartel
% FORMAT spm_dartel_import(job)
% job.matnames  - Names of *_seg_sn.mat files to use
% job.odir      - Output directory
% job.bb        - Bounding box
% job.vox       - Voxel sizes
% job.GM/WM/CSF - Options fo different tissue classes
% job.image     - Options for resliced original image
%
% Rigidly aligned images are generated using info from the seg_sn.mat
% files.  These can be resliced GM, WM or CSF, but also various resliced
% forms of the original image (skull-stripped, bias corrected etc).
%____________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dartel_import.m 5506 2013-05-14 17:13:43Z john $

matnames = job.matnames;
for i=1:numel(matnames),
    p(i) = load(matnames{i});
end;
if numel(p)>0,
    tmp = strvcat(p(1).VG.fname);
    p(1).VG = spm_vol(tmp);
    b0  = spm_load_priors(p(1).VG);
end;
odir = job.odir{1};
bb   = job.bb;
vox  = job.vox;
iopt = job.image;
opt  = [job.GM, job.WM, job.CSF];
for i=1:numel(p),
    preproc_apply(p(i),odir,b0,bb,vox,iopt,opt,matnames{i});
end;

out.cfiles = cell(numel(matnames),numel(opt));
if job.image,
    out.files  = cell(numel(matnames),1);
end
for i=1:numel(matnames),
    [pth,nam,ext] = fileparts(matnames{i});
    nam = nam(1:(numel(nam)-7));
    if job.image,
        fname = fullfile(odir,['r',nam, '.nii']);
        out.files{i} = fname;
    end
    for k1=1:numel(opt),
        if opt(k1),
            fname            = fullfile(odir,['r','c', num2str(k1), nam, '.nii']);
            out.cfiles{i,k1} = fname;
        end
    end
end

return;
%=======================================================================

%=======================================================================
function preproc_apply(p,odir,b0,bb,vx,iopt,opt,matname)
[pth0,nam,ext,num] = spm_fileparts(matname);
[pth ,nam,ext,num] = spm_fileparts(p.VF(1).fname);
P = path_search([nam,ext],{pth,odir,pwd,pth0});
if isempty(P),
    fprintf('Could not find "%s"\n', [nam,ext]);
    return;
end;
p.VF.fname = P;
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

lkp = []; for k=1:Kb, lkp = [lkp ones(1,p.flags.ngaus(k))*k]; end;

spm_progress_bar('init',length(x3),['Working on ' nam],'Planes completed');
M = p.VG(1).mat\p.flags.Affine*p.VF.mat;

if iopt, idat = zeros(d(1:3),'single'); end;
dat = {zeros(d(1:3),'uint8'),zeros(d(1:3),'uint8'),zeros(d(1:3),'uint8')};

for z=1:length(x3),

    % Bias corrected image
    f          = spm_sample_vol(p.VF,x1,x2,o*x3(z),0);
    msk        = (f==0) | ~isfinite(f);
    if ~isempty(bsol),
        cr     = exp(transf(bB1,bB2,bB3(z,:),bsol)).*f;
    else
        cr     = f;
    end

    if iopt,
        if bitand(iopt,2),
            idat(:,:,z) = cr;
        else
            idat(:,:,z) = f;
        end;
    end;

    [t1,t2,t3] = defs(T,z,B1,B2,B3,x1,x2,x3,M);
    q          = zeros([d(1:2) Kb]);
    bg         = ones(d(1:2));
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
    for k1=1:3,
        tmp            = q(:,:,k1);
        tmp            = tmp./sq;
        tmp(msk)       = 0;
        dat{k1}(:,:,z) = uint8(round(255 * tmp));
    end;
    spm_progress_bar('set',z);
end;
spm_progress_bar('clear');

%[dat{1},dat{2},dat{3}] = clean_gwc(dat{1},dat{2},dat{3}, 2); 
if iopt,
    if bitand(iopt,2),
        nwm = 0;
        swm = 0;
        for z=1:numel(x3),
            nwm = nwm + sum(sum(double(dat{2}(:,:,z))));
            swm = swm + sum(sum(double(dat{2}(:,:,z)).*idat(:,:,z)));
        end;
        idat = idat*(double(nwm)/double(swm));
    end;
    if bitand(iopt,4),
        for z=1:numel(x3),
           %msk         = double(dat{1}(:,:,z)) ...
           %            + double(dat{2}(:,:,z)) ...
           %        + 0.5*double(dat{3}(:,:,z));
           %msk         = msk<128;
           %tmp         = idat(:,:,z);
           %tmp(msk)    = 0;
           %idat(:,:,z) = tmp;
           wt           = (double(dat{1}(:,:,z))+double(dat{2}(:,:,z)))/255;
           idat(:,:,z)  = idat(:,:,z).*wt;
        end;
    end;
    for z=1:numel(x3),
        tmp               = idat(:,:,z);
        tmp(~isfinite(double(tmp))) = 0;
        idat(:,:,z)       = tmp;
    end;
end;

% Sort out bounding box etc
[bb1,vx1] = spm_get_bbox(p.VG(1), 'old');
bb(~isfinite(bb)) = bb1(~isfinite(bb));
if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end;
bb(1,:) = vx*ceil(bb(1,:)/vx);
bb(2,:) = vx*floor(bb(2,:)/vx);

% Figure out the mapping from the volumes to create to the original
mm = [
    bb(1,1) bb(1,2) bb(1,3)
    bb(2,1) bb(1,2) bb(1,3)
    bb(1,1) bb(2,2) bb(1,3)    
    bb(2,1) bb(2,2) bb(1,3)
    bb(1,1) bb(1,2) bb(2,3)    
    bb(2,1) bb(1,2) bb(2,3)    
    bb(1,1) bb(2,2) bb(2,3)    
    bb(2,1) bb(2,2) bb(2,3)]';

vx2 = inv(p.VG(1).mat)*[mm ; ones(1,8)];

odim   = abs(round((bb(2,:)-bb(1,:))/vx))+1;
%dimstr = sprintf('%dx%dx%d',odim);
dimstr  = '';

vx1 = [
    1       1       1
    odim(1) 1       1
    1       odim(2) 1
    odim(1) odim(2) 1
    1       1       odim(3)
    odim(1) 1       odim(3)
    1       odim(2) odim(3)
    odim(1) odim(2) odim(3)]';
M      = p.Affine*vx2/[vx1 ; ones(1,8)];
mat0   = p.VF.mat*M;
mat    = [mm ; ones(1,8)]/[vx1 ; ones(1,8)];

fwhm = max(vx./sqrt(sum(p.VF.mat(1:3,1:3).^2))-1,0.01);

for k1=1:numel(opt),
    if opt(k1),
        dat{k1} = decimate(dat{k1},fwhm);
        VT      = struct('fname',fullfile(odir,['r',dimstr,'c', num2str(k1), nam, '.nii']),...
            'dim',  odim,...
            'dt',   [spm_type('uint8') spm_platform('bigend')],...
            'pinfo',[1/255 0]',...
            'mat',mat);
        VT = spm_create_vol(VT);

        Ni             = nifti(VT.fname);
        Ni.mat0        = mat0;
        Ni.mat_intent  = 'Aligned';
        Ni.mat0_intent = 'Aligned';
        create(Ni);

        for i=1:odim(3),
            tmp = spm_slice_vol(dat{k1},M*spm_matrix([0 0 i]),odim(1:2),1)/255;
            VT  = spm_write_plane(VT,tmp,i);
        end;
    end;
end;

if iopt,
        %idat = decimate(idat,fwhm);
        VT      = struct('fname',fullfile(odir,['r',dimstr,nam, '.nii']),...
            'dim',  odim,...
            'dt',   [spm_type('float32') spm_platform('bigend')],...
            'pinfo',[1 0]',...
            'mat',mat);
        VT.descrip = 'Resliced';
        if bitand(iopt,2), VT.descrip = [VT.descrip ', bias corrected']; end;
        if bitand(iopt,4), VT.descrip = [VT.descrip ', skull stripped']; end;
        
        VT = spm_create_vol(VT);

        Ni             = nifti(VT.fname);
        Ni.mat0        = mat0;
        Ni.mat_intent  = 'Aligned';
        Ni.mat0_intent = 'Aligned';
        create(Ni);

        for i=1:odim(3),
            tmp = spm_slice_vol(idat,M*spm_matrix([0 0 i]),odim(1:2),1);
            VT  = spm_write_plane(VT,tmp,i);
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
d2 = [size(T) 1];
t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
t  = B1*t1*B2';
return;
%=======================================================================

%=======================================================================
function dat = decimate(dat,fwhm)
% Convolve the volume in memory (fwhm in voxels).
lim = ceil(2*fwhm);
s  = fwhm/sqrt(8*log(2));
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
function pthnam = path_search(nam,pth)
pthnam = '';
for i=1:numel(pth),
    if exist(fullfile(pth{i},nam),'file'),
        pthnam = fullfile(pth{i},nam);
        return;
    end;
end;

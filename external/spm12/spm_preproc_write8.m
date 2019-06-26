function [cls,M1] = spm_preproc_write8(res,tc,bf,df,mrf,cleanup,bb,vx)
% Write out VBM preprocessed data
% FORMAT [cls,M1] = spm_preproc_write8(res,tc,bf,df,mrf,cleanup,bb,vx)
%__________________________________________________________________________
% Copyright (C) 2008-2016 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_preproc_write8.m 6881 2016-09-19 09:48:54Z john $

% Prior adjustment factor.
% This is a fudge factor to weaken the effects of the tissue priors.  The
% idea here is that the bias from the tissue priors probably needs to be
% reduced because of the spatial smoothing typically used in VBM studies.
% Having the optimal bias/variance tradeoff for each voxel is not the same
% as having the optimal tradeoff for weighted averages over several voxels.

if isfield(res,'mg')
    lkp = res.lkp;
    Kb  = max(lkp);
else
    Kb  = size(res.intensity(1).lik,2);
end
N   = numel(res.image);

if nargin<2, tc = true(Kb,4); end % native, import, warped, warped-mod
if nargin<3, bf = false(N,2); end % field, corrected
if nargin<4, df = false(1,2); end % inverse, forward
if nargin<5, mrf = 1;         end % MRF parameter
if nargin<6, cleanup = 1;     end % Run the ad hoc cleanup

% Read essentials from tpm (it will be cleared later)
tpm = res.tpm;
if ~isstruct(tpm) || ~isfield(tpm, 'bg1'),
    tpm = spm_load_priors8(tpm);
end
d1        = size(tpm.dat{1});
d1        = d1(1:3);
M1        = tpm.M;

% Define orientation and field of view of any "normalised" space
% data that may be generated (wc*.nii, mwc*.nii, rc*.nii & y_*.nii).
if nargin>=7 && any(isfinite(bb(:)))
    % If a bounding box is supplied, combine this with the closest
    % bounding box derived from the dimensions and orientations of
    % the tissue priors.
    if nargin<7, bb = NaN(2,3);   end % Default to TPM bounding box
    if nargin<8, vx = NaN;        end % Default to TPM voxel size
    [bb1,vx1] = spm_get_bbox(tpm.V(1), 'old');
    bb(~isfinite(bb)) = bb1(~isfinite(bb));
    if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end
    bb(1,:) = vx*round(bb(1,:)/vx);
    bb(2,:) = vx*round(bb(2,:)/vx);
    odim    = abs(round((bb(2,1:3)-bb(1,1:3))/vx))+1;

    mm  = [[bb(1,1) bb(1,2) bb(1,3)
            bb(2,1) bb(1,2) bb(1,3)
            bb(1,1) bb(2,2) bb(1,3)
            bb(2,1) bb(2,2) bb(1,3)
            bb(1,1) bb(1,2) bb(2,3)
            bb(2,1) bb(1,2) bb(2,3)
            bb(1,1) bb(2,2) bb(2,3)
            bb(2,1) bb(2,2) bb(2,3)]'; ones(1,8)];
    vx3 = [[1       1       1
            odim(1) 1       1
            1       odim(2) 1
            odim(1) odim(2) 1
            1       1       odim(3)
            odim(1) 1       odim(3)
            1       odim(2) odim(3)
            odim(1) odim(2) odim(3)]'; ones(1,8)];
    mat    = mm/vx3;
else
    % Use the actual dimensions and orientations of
    % the tissue priors.
    odim = tpm.V(1).dim;
    mat  = tpm.V(1).mat;
end


[pth,nam] = fileparts(res.image(1).fname);
ind  = res.image(1).n;
d    = res.image(1).dim(1:3);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3 = 1:d(3);

chan(N) = struct('B1',[],'B2',[],'B3',[],'T',[],'Nc',[],'Nf',[],'ind',[]);
for n=1:N
    d3         = [size(res.Tbias{n}) 1];
    chan(n).B3 = spm_dctmtx(d(3),d3(3),x3);
    chan(n).B2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
    chan(n).B1 = spm_dctmtx(d(1),d3(1),x1(:,1));
    chan(n).T  = res.Tbias{n};

    % Need to fix writing of bias fields or bias corrected images, when the data used are 4D.
    [pth1,nam1] = fileparts(res.image(n).fname);
    chan(n).ind = res.image(n).n;

    if bf(n,2)
        chan(n).Nc      = nifti;
        chan(n).Nc.dat  = file_array(fullfile(pth1,['m', nam1, '.nii']),...
                                     res.image(n).dim(1:3),...
                                     [spm_type('float32') spm_platform('bigend')],...
                                     0,1,0);
        chan(n).Nc.mat  = res.image(n).mat;
        chan(n).Nc.mat0 = res.image(n).mat;
        chan(n).Nc.descrip = 'Bias corrected';
        create(chan(n).Nc);
    end

    if bf(n,1),
        chan(n).Nf      = nifti;
        chan(n).Nf.dat  = file_array(fullfile(pth1,['BiasField_', nam1, '.nii']),...
                                     res.image(n).dim(1:3),...
                                     [spm_type('float32') spm_platform('bigend')],...
                                     0,1,0);
        chan(n).Nf.mat  = res.image(n).mat;
        chan(n).Nf.mat0 = res.image(n).mat;
        chan(n).Nf.descrip = 'Estimated Bias Field';
        create(chan(n).Nf);
    end
end

do_cls   = any(tc(:)) || nargout>=1;
tiss(Kb) = struct('Nt',[]);
for k1=1:Kb
    if tc(k1,4) || any(tc(:,3)) || tc(k1,2) || nargout>=1,
        do_cls  = true;
    end
    if tc(k1,1),
        tiss(k1).Nt      = nifti;
        tiss(k1).Nt.dat  = file_array(fullfile(pth,['c', num2str(k1), nam, '.nii']),...
                                      res.image(n).dim(1:3),...
                                      [spm_type('uint8') spm_platform('bigend')],...
                                      0,1/255,0);
        tiss(k1).Nt.mat  = res.image(n).mat;
        tiss(k1).Nt.mat0 = res.image(n).mat;
        tiss(k1).Nt.descrip = ['Tissue class ' num2str(k1)];
        create(tiss(k1).Nt);
        do_cls = true;
    end
end

prm     = [3 3 3 0 0 0];
Coef    = cell(1,3);
Coef{1} = spm_bsplinc(res.Twarp(:,:,:,1),prm);
Coef{2} = spm_bsplinc(res.Twarp(:,:,:,2),prm);
Coef{3} = spm_bsplinc(res.Twarp(:,:,:,3),prm);

do_defs = any(df);
do_defs = do_cls | do_defs;
if do_defs
    if df(1)
        [pth,nam] =fileparts(res.image(1).fname);
        Ndef      = nifti;
        Ndef.dat  = file_array(fullfile(pth,['iy_', nam, '.nii']),...
                               [res.image(1).dim(1:3),1,3],...
                               [spm_type('float32') spm_platform('bigend')],...
                               0,1,0);
        Ndef.mat  = res.image(1).mat;
        Ndef.mat0 = res.image(1).mat;
        Ndef.descrip = 'Inverse Deformation';
        create(Ndef);
    end
    if df(2) || any(any(tc(:,[2,3,4]))) || nargout>=1,
        y = zeros([res.image(1).dim(1:3),3],'single');
    end
end

spm_progress_bar('init',length(x3),['Working on ' nam],'Planes completed');
M = M1\res.Affine*res.image(1).mat;

if do_cls
    Q = zeros([d(1:3),Kb],'single');
end

for z=1:length(x3)

    % Bias corrected image
    cr = cell(1,N);
    for n=1:N
        f          = spm_sample_vol(res.image(n),x1,x2,o*x3(z),0);
        bf         = exp(transf(chan(n).B1,chan(n).B2,chan(n).B3(z,:),chan(n).T));
        cr{n}      = bf.*f;
        if ~isempty(chan(n).Nc),
            % Write a plane of bias corrected data
            chan(n).Nc.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = cr{n};
        end
        if ~isempty(chan(n).Nf),
            % Write a plane of bias field
            chan(n).Nf.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = bf;
        end
    end


    if do_defs
        % Compute the deformation (mapping voxels in image to voxels in TPM)
        [t1,t2,t3] = defs(Coef,z,res.MT,prm,x1,x2,x3,M);

        if exist('Ndef','var')
            % Write out the deformation to file, adjusting it so mapping is
            % to voxels (voxels in image to mm in TPM)
            tmp = M1(1,1)*t1 + M1(1,2)*t2 + M1(1,3)*t3 + M1(1,4);
            Ndef.dat(:,:,z,1,1) = tmp;
            tmp = M1(2,1)*t1 + M1(2,2)*t2 + M1(2,3)*t3 + M1(2,4);
            Ndef.dat(:,:,z,1,2) = tmp;
            tmp = M1(3,1)*t1 + M1(3,2)*t2 + M1(3,3)*t3 + M1(3,4);
            Ndef.dat(:,:,z,1,3) = tmp;
        end

        if exist('y','var')
            % If needed later, save in variable y
            y(:,:,z,1) = t1;
            y(:,:,z,2) = t2;
            y(:,:,z,3) = t3;
        end

        if do_cls
            % Generate variable Q if tissue classes are needed
           %msk = (f==0) | ~isfinite(f);

            if isfield(res,'mg')
                % Parametric representation of intensity distributions
                q   = zeros([d(1:2) Kb]);
                q1  = likelihoods(cr,[],res.mg,res.mn,res.vr);
                q1  = reshape(q1,[d(1:2),numel(res.mg)]);
                b   = spm_sample_priors8(tpm,t1,t2,t3);
                wp  = res.wp;
                s   = zeros(size(b{1}));
                for k1 = 1:Kb
                    b{k1} = wp(k1)*b{k1};
                    s     = s + b{k1};
                end
                for k1=1:Kb
                    q(:,:,k1) = sum(q1(:,:,lkp==k1),3).*(b{k1}./s);
                end
            else
                % Nonparametric representation of intensity distributions
                q   = spm_sample_priors8(tpm,t1,t2,t3);
                wp  = res.wp;
                s   = zeros(size(q{1}));
                for k1 = 1:Kb
                    q{k1} = wp(k1)*q{k1};
                    s     = s + q{k1};
                end
                for k1 = 1:Kb
                    q{k1} = q{k1}./s;
                end
                q   = cat(3,q{:});

                for n=1:N
                    tmp = round(cr{n}*res.intensity(n).interscal(2) + res.intensity(n).interscal(1));
                    tmp = min(max(tmp,1),size(res.intensity(n).lik,1));
                    for k1=1:Kb
                        likelihood = res.intensity(n).lik(:,k1);
                        q(:,:,k1)  = q(:,:,k1).*likelihood(tmp);
                    end
                end
            end
            Q(:,:,z,:) = reshape(q,[d(1:2),1,Kb]);
        end
    end
    spm_progress_bar('set',z);
end
spm_progress_bar('clear');


cls   = cell(1,Kb);
if do_cls
    P = zeros([d(1:3),Kb],'uint8');
    if mrf==0
        % Normalise to sum to 1
        sQ = (sum(Q,4)+eps)/255;
        for k1=1:size(Q,4)
            P(:,:,:,k1) = uint8(round(Q(:,:,:,k1)./sQ));
        end
        clear sQ
    else
        % Use a MRF cleanup procedure
        nmrf_its = 10;
        spm_progress_bar('init',nmrf_its,['MRF: Working on ' nam],'Iterations completed');
        G   = ones([Kb,1],'single')*mrf;
        vx2 = 1./single(sum(res.image(1).mat(1:3,1:3).^2));
        %save PQG P Q G tiss Kb x3 ind
        for iter=1:nmrf_its
            spm_mrf(P,Q,G,vx2);
            spm_progress_bar('set',iter);
        end
    end
    clear Q

    if cleanup
        % Use an ad hoc brain cleanup procedure
        if size(P,4)>5
            P = clean_gwc(P,cleanup);
        else
            warning('Cleanup not done.');
        end
    end

    % Write tissues if necessary
    for k1=1:Kb
        if ~isempty(tiss(k1).Nt)
            for z=1:length(x3)
                tmp = double(P(:,:,z,k1))/255;
                tiss(k1).Nt.dat(:,:,z,ind(1),ind(2)) = tmp;
            end
        end
    end
    spm_progress_bar('clear');

    % Put tissue classes into a cell array...
    for k1=1:Kb
        if tc(k1,4) || any(tc(:,3)) || tc(k1,2) || nargout>=1
            cls{k1} = P(:,:,:,k1);
        end
    end
    clear P % ...and remove the original 4D array
end

clear tpm
M0  = res.image(1).mat;

if any(tc(:,2))
    % "Imported" tissue class images

    % Generate mm coordinates of where deformations map from
    x      = affind(rgrid(d),M0);

    % Generate mm coordinates of where deformation maps to
    y1     = affind(y,M1);

    % Procrustes analysis to compute the closest rigid-body
    % transformation to the deformation, weighted by the
    % interesting tissue classes.
    ind    = find(tc(:,2)); % Saved tissue classes
    [dummy,R]  = spm_get_closest_affine(x,y1,single(cls{ind(1)})/255);
    clear x y1

    mat0   = R\mat; % Voxel-to-world of original image space

    [pth,nam] = fileparts(res.image(1).fname);
    fwhm   = max(vx./sqrt(sum(res.image(1).mat(1:3,1:3).^2))-1,0.01);
    for k1=1:size(tc,1)
        if tc(k1,2)

            % Low pass filtering to reduce aliasing effects in downsampled images,
            % then reslice and write to disk
            tmp1    = decimate(single(cls{k1}),fwhm);
            Ni      = nifti;
            Ni.dat  = file_array(fullfile(pth,['rc', num2str(k1), nam, '.nii']),...
                                 odim,...
                                 [spm_type('float32') spm_platform('bigend')],...
                                 0,1,0);
            Ni.mat         = mat;
            Ni.mat_intent  = 'Aligned';
            Ni.mat0        = mat0;
            Ni.mat0_intent = 'Aligned';
            Ni.descrip     = ['Imported Tissue ' num2str(k1)];
            create(Ni);

            for i=1:odim(3)
                tmp           = spm_slice_vol(tmp1,M0\mat0*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255;
                Ni.dat(:,:,i) = tmp;
            end
            clear tmp1
        end
    end
end

if any(tc(:,3)) || any(tc(:,4)) || nargout>=1 || df(2)
    % Adjust stuff so that warped data (and deformations) have the
    % desired bounding box and voxel sizes, instead of being the same
    % as those of the tissue probability maps.
    M = mat\M1;
    for i=1:size(y,3)
        t1         = y(:,:,i,1);
        t2         = y(:,:,i,2);
        t3         = y(:,:,i,3);
        y(:,:,i,1) = M(1,1)*t1 + M(1,2)*t2 + M(1,3)*t3 + M(1,4);
        y(:,:,i,2) = M(2,1)*t1 + M(2,2)*t2 + M(2,3)*t3 + M(2,4);
        y(:,:,i,3) = M(3,1)*t1 + M(3,2)*t2 + M(3,3)*t3 + M(3,4);
    end
    M1 = mat;
    d1 = odim;
end


if any(tc(:,3)) || any(tc(:,4)) || nargout>=1

    if any(tc(:,3))
        C = zeros([d1,Kb],'single');
    end

    spm_progress_bar('init',Kb,'Warped Tissue Classes','Classes completed');
    for k1 = 1:Kb
        if ~isempty(cls{k1})
            c = single(cls{k1})/255;
            if any(tc(:,3))
                [c,w]       = spm_diffeo('push',c,y,d1(1:3));
                vx          = sqrt(sum(M1(1:3,1:3).^2));
                spm_field('boundary',1);
                C(:,:,:,k1) = spm_field(w,c,[vx  1e-6 1e-4 0  3 2]);
                clear w
            else
                c      = spm_diffeo('push',c,y,d1(1:3));
            end
            if nargout>=1
                cls{k1} = c;
            end
            if tc(k1,4)
                N      = nifti;
                N.dat  = file_array(fullfile(pth,['mwc', num2str(k1), nam, '.nii']),...
                                    d1,...
                                    [spm_type('float32') spm_platform('bigend')],...
                                    0,1,0);
                N.mat  = M1;
                N.mat0 = M1;
                N.descrip = ['Jac. sc. warped tissue class ' num2str(k1)];
                create(N);
                N.dat(:,:,:) = c*abs(det(M0(1:3,1:3))/det(M1(1:3,1:3)));
            end
            spm_progress_bar('set',k1);
        end
    end
    spm_progress_bar('Clear');

    if any(tc(:,3))
        spm_progress_bar('init',Kb,'Writing Warped Tis Cls','Classes completed');
        C = max(C,eps);
        s = sum(C,4);
        for k1=1:Kb
            if tc(k1,3)
                N      = nifti;
                N.dat  = file_array(fullfile(pth,['wc', num2str(k1), nam, '.nii']),...
                                    d1,'uint8',0,1/255,0);
                N.mat  = M1;
                N.mat0 = M1;
                N.descrip = ['Warped tissue class ' num2str(k1)];
                create(N);
                N.dat(:,:,:) = C(:,:,:,k1)./s;
            end
            spm_progress_bar('set',k1);
        end
        spm_progress_bar('Clear');
        clear C s
    end
end

if df(2)
    y         = spm_diffeo('invdef',y,d1,eye(4),M0);
    y         = spm_extrapolate_def(y,M1);
    N         = nifti;
    N.dat     = file_array(fullfile(pth,['y_', nam, '.nii']),...
                           [d1,1,3],'float32',0,1,0);
    N.mat     = M1;
    N.mat0    = M1;
    N.descrip = 'Deformation';
    create(N);
    N.dat(:,:,:,:,:) = reshape(y,[d1,1,3]);
end

return;


%==========================================================================
% function [x1,y1,z1] = defs(sol,z,MT,prm,x0,y0,z0,M)
%==========================================================================
function [x1,y1,z1] = defs(sol,z,MT,prm,x0,y0,z0,M)
iMT = inv(MT);
x1  = x0*iMT(1,1)+iMT(1,4);
y1  = y0*iMT(2,2)+iMT(2,4);
z1  = (z0(z)*iMT(3,3)+iMT(3,4))*ones(size(x1));
x1a = x0    + spm_bsplins(sol{1},x1,y1,z1,prm);
y1a = y0    + spm_bsplins(sol{2},x1,y1,z1,prm);
z1a = z0(z) + spm_bsplins(sol{3},x1,y1,z1,prm);
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
return;


%==========================================================================
% function t = transf(B1,B2,B3,T)
%==========================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t = zeros(size(B1,1),size(B2,1),size(B3,1));
end
return;


%==========================================================================
% function p = likelihoods(f,bf,mg,mn,vr)
%==========================================================================
function p = likelihoods(f,bf,mg,mn,vr)
K  = numel(mg);
N  = numel(f);
M  = numel(f{1});
cr = zeros(M,N);
for n=1:N
    if isempty(bf)
        cr(:,n) = double(f{n}(:));
    else
        cr(:,n) = double(f{n}(:).*bf{n}(:));
    end
end
p  = ones(numel(f{1}),K);
for k=1:K
    amp    = mg(k)/sqrt((2*pi)^N * det(vr(:,:,k)));
    d      = bsxfun(@minus,cr,mn(:,k)')/chol(vr(:,:,k));
    p(:,k) = amp*exp(-0.5*sum(d.*d,2)) + eps;
end
return;


%==========================================================================
% function dat = decimate(dat,fwhm)
%==========================================================================
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


%==========================================================================
% function y1 = affind(y0,M)
%==========================================================================
function y1 = affind(y0,M)
y1 = zeros(size(y0),'single');
for d=1:3
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1) + y0(:,:,:,2)*M(d,2) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
return;


%==========================================================================
% function x = rgrid(d)
%==========================================================================
function x = rgrid(d)
x = zeros([d(1:3) 3],'single');
[x1,x2] = ndgrid(single(1:d(1)),single(1:d(2)));
for i=1:d(3)
    x(:,:,i,1) = x1;
    x(:,:,i,2) = x2;
    x(:,:,i,3) = single(i);
end
return;


%==========================================================================
% function [P] = clean_gwc(P,level)
%==========================================================================
function [P] = clean_gwc(P,level)
if nargin<2, level = 1; end

b    = P(:,:,:,2);

% Build a 3x3x3 seperable smoothing kernel
%--------------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

th1 = 0.15;
if level==2, th1 = 0.2; end
% Erosions and conditional dilations
%--------------------------------------------------------------------------
niter  = 32;
niter2 = 32;
spm_progress_bar('Init',niter+niter2,'Extracting Brain','Iterations completed');
for j=1:niter
    if j>2, th=th1; else th=0.6; end  % Dilate after two its of erosion
    for i=1:size(b,3)
        gp       = double(P(:,:,i,1));
        wp       = double(P(:,:,i,2));
        bp       = double(b(:,:,i))/255;
        bp       = (bp>th).*(wp+gp);
        b(:,:,i) = uint8(round(bp));
    end
    spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
    spm_progress_bar('Set',j);
end

% Also clean up the CSF.
if niter2 > 0
    c = b;
    for j=1:niter2
        for i=1:size(b,3)
            gp       = double(P(:,:,i,1));
            wp       = double(P(:,:,i,2));
            cp       = double(P(:,:,i,3));
            bp       = double(c(:,:,i))/255;
            bp       = (bp>th).*(wp+gp+cp);
            c(:,:,i) = uint8(round(bp));
        end
        spm_conv_vol(c,c,kx,ky,kz,-[1 1 1]);
        spm_progress_bar('Set',j+niter);
    end
end

th = 0.05;
for i=1:size(b,3)
    slices = cell(1,size(P,4));
    for k1=1:size(P,4)
        slices{k1} = double(P(:,:,i,k1))/255;
    end
    bp        = double(b(:,:,i))/255;
    bp        = ((bp>th).*(slices{1}+slices{2}))>th;
    slices{1} = slices{1}.*bp;
    slices{2} = slices{2}.*bp;

    if niter2>0
        cp        = double(c(:,:,i))/255;
        cp        = ((cp>th).*(slices{1}+slices{2}+slices{3}))>th;
        slices{3} = slices{3}.*cp;
    end
    slices{5} = slices{5}+1e-4; % Add a little to the soft tissue class
    tot       = zeros(size(bp))+eps;
    for k1=1:size(P,4)
        tot   = tot + slices{k1};
    end
    for k1=1:size(P,4)
        P(:,:,i,k1) = uint8(round(slices{k1}./tot*255));
    end 
end
spm_progress_bar('Clear');

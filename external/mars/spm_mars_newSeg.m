function cls = spm_mars_newSeg(res,tc,bf,df,mode,tcm)
% FORMAT cls = spm_mars_newSeg(res,tc,bf,df,mode,tcm)
%
% Write out the results from toolbox 'New Segment'. This function is only for
% the purpose of making toolbox 'MARS' backwards compatible with 'New Segment'.
%
% Adapted from spm_preproc_write8 by John Ashburner
% $Id: spm_preproc_write8.m 4337 2011-05-31 16:59:44Z john $
% Copyright (C) 2008 Wellcome Department of Imaging Neuroscience
%
% Modified by Yu (Andy) Huang
% $Id: spm_mars_newSeg.m 2015-07-27 andy$
% Neural Engineering Lab, Dept. of Biomedical Engineering, City College of New York
% yhuang16@citymail.cuny.edu

% Read essentials from tpm (it will be cleared later)
tpm = res.tpm;
if ~isstruct(tpm) || ~isfield(tpm, 'bg1'),
    tpm = spm_load_priors8(tpm);
end
d1        = size(tpm.dat{1});
d1        = d1(1:3);
M1        = tpm.M;
[bb1 vx1] = spm_get_bbox(tpm.V(1), 'old');

if isfield(res,'mg'),
    lkp = res.lkp;
    Kb  = max(lkp);
else
    Kb  = size(res.intensity(1).lik,2);
end

N   = numel(res.image);
if nargin<2, tc = true(Kb,4); end % native, import, warped, warped-mod
if nargin<3, bf = true(N,2);  end % field, corrected
if nargin<4, df = true(1,2);  end % inverse, forward
% if nargin<5, mrf= 2;          end % MRF parameter  % andy 2013-05-03

[pth,nam]=fileparts(res.image(1).fname);
ind  = res.image(1).n;
d    = res.image(1).dim(1:3);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3  = 1:d(3);

chan(N) = struct('B1',[],'B2',[],'B3',[],'T',[],'Nc',[],'Nf',[],'ind',[]);
for n=1:N,
    d3         = [size(res.Tbias{n}) 1];
    chan(n).B3 = spm_dctmtx(d(3),d3(3),x3);
    chan(n).B2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
    chan(n).B1 = spm_dctmtx(d(1),d3(1),x1(:,1));
    chan(n).T  = res.Tbias{n};

    [pth1,nam1,ext1] = fileparts(res.image(n).fname);
    chan(n).ind      = res.image(n).n;

    if bf(n,2),
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

do_cls   = any(tc(:)) || nargout>1;
tiss(Kb) = struct('Nt',[]);
for k1=1:Kb,
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
    end;
end

prm     = [3 3 3 0 0 0];
Coef    = cell(1,3);
Coef{1} = spm_bsplinc(res.Twarp(:,:,:,1),prm);
Coef{2} = spm_bsplinc(res.Twarp(:,:,:,2),prm);
Coef{3} = spm_bsplinc(res.Twarp(:,:,:,3),prm);

do_defs = any(df);
do_defs = do_cls | do_defs;
if do_defs,
    if df(1),
        [pth,nam,ext1]=fileparts(res.image(1).fname);
        Ndef      = nifti;
        Ndef.dat  = file_array(fullfile(pth,['iy_', nam1, '.nii']),...
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

for z=1:length(x3),

    % Bias corrected image
    cr = cell(1,N);
    for n=1:N,
        f          = spm_sample_vol(res.image(n),x1,x2,o*x3(z),0);
        bf         = exp(transf(chan(n).B1,chan(n).B2,chan(n).B3(z,:),chan(n).T));
        cr{n}      = bf.*f;
        if ~isempty(chan(n).Nc),
            % Write a plane of bias corrected data
            chan(n).Nc.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = cr{n};
        end;
        if ~isempty(chan(n).Nf),
            % Write a plane of bias field
            chan(n).Nf.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = bf;
        end;
    end


    if do_defs,
        [t1,t2,t3] = defs(Coef,z,res.MT,prm,x1,x2,x3,M);
        if exist('Ndef','var'),
            tmp = M1(1,1)*t1 + M1(1,2)*t2 + M1(1,3)*t3 + M1(1,4);
            Ndef.dat(:,:,z,1,1) = tmp;
            tmp = M1(2,1)*t1 + M1(2,2)*t2 + M1(2,3)*t3 + M1(2,4);
            Ndef.dat(:,:,z,1,2) = tmp;
            tmp = M1(3,1)*t1 + M1(3,2)*t2 + M1(3,3)*t3 + M1(3,4);
            Ndef.dat(:,:,z,1,3) = tmp;
        end

        if exist('y','var'),
            y(:,:,z,1) = t1;
            y(:,:,z,2) = t2;
            y(:,:,z,3) = t3;
        end

        if do_cls,
            msk = (f==0) | ~isfinite(f);

            if isfield(res,'mg'),
                q   = zeros([d(1:2) Kb]);
                q1  = likelihoods(cr,[],res.mg,res.mn,res.vr);
                q1  = reshape(q1,[d(1:2),numel(res.mg)]);
                b   = spm_sample_priors8(tpm,t1,t2,t3);
                for k1=1:Kb,
                    q(:,:,k1) = sum(q1(:,:,lkp==k1),3).*b{k1};
                end
            else
                q   = spm_sample_priors8(tpm,t1,t2,t3);
                q   = cat(3,q{:});
                for n=1:N,
                    tmp = round(cr{n}*res.intensity(n).interscal(2) + res.intensity(n).interscal(1));
                    tmp = min(max(tmp,1),size(res.intensity(n).lik,1));
                    for k1=1:Kb,
                        likelihood = res.intensity(n).lik(:,k1);
                        q(:,:,k1)  = q(:,:,k1).*likelihood(tmp);
                    end
                end
            end
            Q(:,:,z,:) = reshape(q,[d(1:2),1,Kb]); % final probabilities not normalized

        end
    end
    spm_progress_bar('set',z);
end
spm_progress_bar('clear');

cls   = cell(1,Kb);
if do_cls
    P = zeros([d(1:3),Kb],'uint8');
    for z=1:length(x3),
        sq = sum(Q(:,:,z,:),4) + eps^2;
        for k1=1:Kb
            P(:,:,z,k1) = uint8(round(255 * Q(:,:,z,k1)./sq)); % nomalized probabilities
        end
    end
%     if mrf~=0, nmrf_its = 10; else nmrf_its = 1; end % andy 2013-05-03

%     if mrf==0, % Done this way as spm_mrf is not compiled for all platforms yet
    if mode == 0  % andy 2013-05-03
        sQ = (sum(Q,4)+eps)/255;
        for k1=1:size(Q,4)
            P(:,:,:,k1) = uint8(round(Q(:,:,:,k1)./sQ));
        end
        clear sQ
    else
        nmrf_its = 10; % andy 2013-05-03
        spm_progress_bar('init',nmrf_its,['MRF: Working on ' nam],'Iterations completed');
        vx2 = single(sum(res.image(1).mat(1:3,1:3).^2));
        if iscellstr(tcm) % if user specifies TCM, load it  % andy 2013-05-03
            load(tcm{1},'C');
            if ndims(C)~=2
                warning('Warn:convert',...
                    'Local or regional TCM detected. For MRF-based clean-up provided by SPM8, user-provided TCM must be a global TCM!\nSPM8 default config is used instead...\n');
                G   = ones([Kb,1],'single')*2;
            else
            G = single(C);  % andy 2013-04-22
            end
        else
            G   = ones([Kb,1],'single')*2; % andy 2013-05-03
            %save PQG P Q G tiss Kb x3 ind
        end
        
        for iter=1:nmrf_its,
            spm_mrf(P,Q,G,vx2);
            spm_progress_bar('set',iter);
        end
    end

    clear Q

    for k1=1:Kb,
        if ~isempty(tiss(k1).Nt),
            for z=1:length(x3),
                tmp = double(P(:,:,z,k1))/255;
                tiss(k1).Nt.dat(:,:,z,ind(1),ind(2)) = tmp;
            end
        end
    end
    spm_progress_bar('clear');

    for k1=1:Kb,
        if tc(k1,4) || any(tc(:,3)) || tc(k1,2) || nargout>=1,
            cls{k1} = P(:,:,:,k1);
        end
    end
    clear P
end

clear tpm
M0 = res.image(1).mat;

if any(tc(:,2)),

    bb = nan(2,3);
    vx = 1.5;
    % Sort out bounding box etc
    bb(~isfinite(bb)) = bb1(~isfinite(bb));
    if ~isfinite(vx), vx = abs(prod(vx1))^(1/3); end;
    bb(1,:) = vx*round(bb(1,:)/vx);
    bb(2,:) = vx*round(bb(2,:)/vx);

    % Figure out the mapping from the volumes to create to the original
    mm = [[
        bb(1,1) bb(1,2) bb(1,3)
        bb(2,1) bb(1,2) bb(1,3)
        bb(1,1) bb(2,2) bb(1,3)
        bb(2,1) bb(2,2) bb(1,3)
        bb(1,1) bb(1,2) bb(2,3)
        bb(2,1) bb(1,2) bb(2,3)
        bb(1,1) bb(2,2) bb(2,3)
        bb(2,1) bb(2,2) bb(2,3)]'; ones(1,8)];

    vx2  = M1\mm;
    odim = abs(round((bb(2,1:3)-bb(1,1:3))/vx))+1;
    vx3  = [[
        1       1       1
        odim(1) 1       1
        1       odim(2) 1
        odim(1) odim(2) 1
        1       1       odim(3)
        odim(1) 1       odim(3)
        1       odim(2) odim(3)
        odim(1) odim(2) odim(3)]'; ones(1,8)];

    x      = affind(rgrid(d),M0);
    y1     = affind(y,M1);
    ind    = find(tc(:,2));
    [M,R]  = spm_get_closest_affine(x,y1,single(cls{ind(1)})/255);
    clear x y1

    M      = M0\inv(R)*M1*vx2/vx3;
    mat0   =         R\M1*vx2/vx3;
    mat    = mm/vx3;

    fwhm = max(vx./sqrt(sum(res.image(1).mat(1:3,1:3).^2))-1,0.01);
    for k1=1:size(tc,1),
        if tc(k1,2),
            tmp1     = decimate(single(cls{k1}),fwhm);
            [pth,nam,ext1]=fileparts(res.image(1).fname);
            VT      = struct('fname',fullfile(pth,['rc', num2str(k1), nam, '.nii']),...
                'dim',  odim,...
                'dt',   [spm_type('float32') spm_platform('bigend')],...
                'pinfo',[1.0 0]',...
                'mat',mat);
            VT = spm_create_vol(VT);

            Ni             = nifti(VT.fname);
            Ni.mat0        = mat0;
            Ni.mat_intent  = 'Aligned';
            Ni.mat0_intent = 'Aligned';
            create(Ni);

            for i=1:odim(3),
                tmp = spm_slice_vol(tmp1,M*spm_matrix([0 0 i]),odim(1:2),[1,NaN])/255;
                VT  = spm_write_plane(VT,tmp,i);
            end
            clear tmp1
        end
    end
end

if any(tc(:,3)),
    C = zeros([d1,Kb],'single');
end

if any(tc(:,3)) || any(tc(:,4)) || nargout>=1,
    spm_progress_bar('init',Kb,'Warped Tissue Classes','Classes completed');
    for k1 = 1:Kb,
        if ~isempty(cls{k1}),
            c = single(cls{k1})/255;
            if any(tc(:,3)),
                [c,w]  = dartel3('push',c,y,d1(1:3));
                vx          = sqrt(sum(M1(1:3,1:3).^2));
                C(:,:,:,k1) = optimNn(w,c,[1  vx  1e-4 1e-6 0  3 2]);
                clear w
            else
                c      = dartel3('push',c,y,d1(1:3));
            end
            if nargout>=1,
                cls{k1} = c;
            end
            if tc(k1,4),
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
end

if any(tc(:,3)),
    spm_progress_bar('init',Kb,'Writing Warped Tis Cls','Classes completed');
    C = max(C,eps);
    s = sum(C,4);
    for k1=1:Kb,
        if tc(k1,3),
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

if df(2),
    y         = spm_invert_def(y,M1,d1,M0,[1 0]);
    N         = nifti;
    N.dat     = file_array(fullfile(pth,['y_', nam1, '.nii']),...
                           [d1,1,3],'float32',0,1,0);
    N.mat     = M1;
    N.mat0    = M1;
    N.descrip = 'Deformation';
    create(N);
    N.dat(:,:,:,:,:) = reshape(y,[d1,1,3]);
end

return;
%=======================================================================

%=======================================================================
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
function p = likelihoods(f,bf,mg,mn,vr)
K  = numel(mg);
N  = numel(f);
M  = numel(f{1});
cr = zeros(M,N);
for n=1:N,
    if isempty(bf),
        cr(:,n) = double(f{n}(:));
    else
        cr(:,n) = double(f{n}(:).*bf{n}(:));
    end
end
p  = ones(numel(f{1}),K);
for k=1:K,
    amp    = mg(k)/sqrt((2*pi)^N * det(vr(:,:,k)));
    d      = cr - repmat(mn(:,k)',M,1);
    p(:,k) = amp * exp(-0.5* sum(d.*(d/vr(:,:,k)),2));
end
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
function y1 = affind(y0,M)
y1 = zeros(size(y0),'single');
for d=1:3,
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1) + y0(:,:,:,2)*M(d,2) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
%=======================================================================

%=======================================================================
function x = rgrid(d)
x = zeros([d(1:3) 3],'single');
[x1,x2] = ndgrid(single(1:d(1)),single(1:d(2)));
for i=1:d(3),
    x(:,:,i,1) = x1;
    x(:,:,i,2) = x2;
    x(:,:,i,3) = single(i);
end
%=======================================================================


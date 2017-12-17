function out = spm_deformations(job)
% Various deformation field utilities
% FORMAT out = spm_deformations(job)
% job - a job created via spm_cfg_deformations.m
% out - a struct with fields
%       .def    - file name of created deformation field
%       .warped - file names of warped images
%
% See spm_cfg_deformations.m for more information.
%__________________________________________________________________________
% Copyright (C) 2005-2015 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_deformations.m 6577 2015-10-15 15:22:11Z volkmar $


[Def,mat] = get_comp(job.comp);
out = struct('def',{{}},'warped',{{}},'surf',{{}},'jac',{{}});
for i=1:numel(job.out)
    fn = fieldnames(job.out{i});
    fn = fn{1};
    switch fn
    case 'savedef'
        out.def    = [out.def;    save_def(Def,mat,job.out{i}.(fn))];
    case 'pull'
        out.warped = [out.warped; pull_def(Def,mat,job.out{i}.(fn))];
    case 'push'
        out.warped = [out.warped; push_def(Def,mat,job.out{i}.(fn))];
    case 'surf'
        out.surf   = [out.surf;   surf_def(Def,mat,job.out{i}.(fn))];
    case 'savejac'
        out.jac    = [out.jac;    jac_def(Def,mat,job.out{i}.(fn))];
    otherwise
        error('Unknown option');
    end
end


%==========================================================================
% function [Def,mat] = get_comp(job)
%==========================================================================
function [Def,mat] = get_comp(job)
% Return the composition of a number of deformation fields.

if isempty(job)
    error('Empty list of jobs in composition');
end
[Def,mat] = get_job(job{1});
for i=2:numel(job)
    Def1         = Def;
    mat1         = mat;
    [Def,mat]    = get_job(job{i});
    M            = inv(mat1);
    tmp          = zeros(size(Def),'single');
    tmp(:,:,:,1) = M(1,1)*Def(:,:,:,1)+M(1,2)*Def(:,:,:,2)+M(1,3)*Def(:,:,:,3)+M(1,4);
    tmp(:,:,:,2) = M(2,1)*Def(:,:,:,1)+M(2,2)*Def(:,:,:,2)+M(2,3)*Def(:,:,:,3)+M(2,4);
    tmp(:,:,:,3) = M(3,1)*Def(:,:,:,1)+M(3,2)*Def(:,:,:,2)+M(3,3)*Def(:,:,:,3)+M(3,4);
    Def(:,:,:,1) = single(spm_diffeo('bsplins',Def1(:,:,:,1),tmp,[1,1,1,0,0,0]));
    Def(:,:,:,2) = single(spm_diffeo('bsplins',Def1(:,:,:,2),tmp,[1,1,1,0,0,0]));
    Def(:,:,:,3) = single(spm_diffeo('bsplins',Def1(:,:,:,3),tmp,[1,1,1,0,0,0]));
    clear tmp
end


%==========================================================================
% function [Def,mat] = get_job(job)
%==========================================================================
function [Def,mat] = get_job(job)
% Determine what is required, and pass the relevant bit of the
% job out to the appropriate function.

fn = fieldnames(job);
fn = fn{1};
switch fn
    case {'comp'}
        [Def,mat] = get_comp(job.(fn));
    case {'def'}
        [Def,mat] = get_def(job.(fn));
    case {'dartel'}
        [Def,mat] = get_dartel(job.(fn));
    case {'sn2def'}
        [Def,mat] = get_sn2def(job.(fn));
    case {'inv'}
        [Def,mat] = get_inv(job.(fn));
    case {'id'}
        [Def,mat] = get_id(job.(fn));
    case {'idbbvox'}
        [Def,mat] = get_idbbvox(job.(fn));
    otherwise
        error('Unrecognised job type');
end


%==========================================================================
% function [Def,mat] = get_sn2def(job)
%==========================================================================
function [Def,mat] = get_sn2def(job)
% Convert a SPM _sn.mat file into a deformation field, and return it.

vox = job.vox;
bb  = job.bb;
sn  = load(job.matname{1});

if any(isfinite(bb(:))) || any(isfinite(vox))
    [bb0, vox0] = spm_get_bbox(sn.VG(1));
    if any(~isfinite(vox)), vox = vox0; end
    if any(~isfinite(bb)),  bb  = bb0;  end
    bb  = sort(bb);
    vox = abs(vox);

    % Adjust bounding box slightly - so it rounds to closest voxel.
    bb(:,1) = round(bb(:,1)/vox(1))*vox(1);
    bb(:,2) = round(bb(:,2)/vox(2))*vox(2);
    bb(:,3) = round(bb(:,3)/vox(3))*vox(3);

    M   = sn.VG(1).mat;
    vxg = sqrt(sum(M(1:3,1:3).^2));
    ogn = M\[0 0 0 1]';
    ogn = ogn(1:3)';

    % Convert range into range of voxels within template image
    x   = (bb(1,1):vox(1):bb(2,1))/vxg(1) + ogn(1);
    y   = (bb(1,2):vox(2):bb(2,2))/vxg(2) + ogn(2);
    z   = (bb(1,3):vox(3):bb(2,3))/vxg(3) + ogn(3);

    og  = -vxg.*ogn;
    of  = -vox.*(round(-bb(1,:)./vox)+1);
    M1  = [vxg(1) 0 0 og(1) ; 0 vxg(2) 0 og(2) ; 0 0 vxg(3) og(3) ; 0 0 0 1];
    M2  = [vox(1) 0 0 of(1) ; 0 vox(2) 0 of(2) ; 0 0 vox(3) of(3) ; 0 0 0 1];
    mat = sn.VG(1).mat*(M1\M2);
    % dim = [length(x) length(y) length(z)];
else
    dim = sn.VG(1).dim;
    x   = 1:dim(1);
    y   = 1:dim(2);
    z   = 1:dim(3);
    mat = sn.VG(1).mat;
end

[X,Y] = ndgrid(x,y);

st = size(sn.Tr);

if (prod(st) == 0)
    affine_only = true;
    basX = 0;
    basY = 0;
    basZ = 0;
else
    affine_only = false;
    basX = spm_dctmtx(sn.VG(1).dim(1),st(1),x-1);
    basY = spm_dctmtx(sn.VG(1).dim(2),st(2),y-1);
    basZ = spm_dctmtx(sn.VG(1).dim(3),st(3),z-1);
end

Def = zeros([numel(x),numel(y),numel(z),3],'single');

for j=1:length(z)
    if (~affine_only)
        tx = reshape( reshape(sn.Tr(:,:,:,1),st(1)*st(2),st(3)) *basZ(j,:)', st(1), st(2) );
        ty = reshape( reshape(sn.Tr(:,:,:,2),st(1)*st(2),st(3)) *basZ(j,:)', st(1), st(2) );
        tz = reshape( reshape(sn.Tr(:,:,:,3),st(1)*st(2),st(3)) *basZ(j,:)', st(1), st(2) );

        X1 = X    + basX*tx*basY';
        Y1 = Y    + basX*ty*basY';
        Z1 = z(j) + basX*tz*basY';
    end

    Mult = sn.VF.mat*sn.Affine;
    if (~affine_only)
        X2= Mult(1,1)*X1 + Mult(1,2)*Y1 + Mult(1,3)*Z1 + Mult(1,4);
        Y2= Mult(2,1)*X1 + Mult(2,2)*Y1 + Mult(2,3)*Z1 + Mult(2,4);
        Z2= Mult(3,1)*X1 + Mult(3,2)*Y1 + Mult(3,3)*Z1 + Mult(3,4);
    else
        X2= Mult(1,1)*X + Mult(1,2)*Y + (Mult(1,3)*z(j) + Mult(1,4));
        Y2= Mult(2,1)*X + Mult(2,2)*Y + (Mult(2,3)*z(j) + Mult(2,4));
        Z2= Mult(3,1)*X + Mult(3,2)*Y + (Mult(3,3)*z(j) + Mult(3,4));
    end

    Def(:,:,j,1) = single(X2);
    Def(:,:,j,2) = single(Y2);
    Def(:,:,j,3) = single(Z2);
end


%==========================================================================
% function [Def,mat] = get_def(job)
%==========================================================================
function [Def,mat] = get_def(job)
% Load a deformation field saved as an image
Nii = nifti(job{1});
Def = single(Nii.dat(:,:,:,1,:));
d   = size(Def);
if d(4)~=1 || d(5)~=3, error('Deformation field is wrong!'); end
Def = reshape(Def,[d(1:3) d(5)]);
mat = Nii.mat;

%==========================================================================
% function [Def,mat] = get_dartel(job)
%==========================================================================
function [Def,mat] = get_dartel(job)

Nii    = nifti(job.flowfield{1});

if ~isempty(job.template{1})
   %Nt        = nifti(job.template{1});
    [pth,nam] = fileparts(job.template{1});
    if exist(fullfile(pth,[nam '_2mni.mat']),'file')
        load(fullfile(pth,[nam '_2mni.mat']),'mni');
    else
        % Affine registration of Dartel Template with MNI space.
        %------------------------------------------------------------------
        fprintf('** Affine registering "%s" with MNI space **\n', nam);
        tpm        = fullfile(spm('Dir'),'tpm','TPM.nii');
        Mmni       = spm_get_space(tpm);
        Nt         = nifti(job.template{1});
        mni.affine = Mmni/spm_klaff(Nt,tpm);
        mni.code   = 'MNI152';
        save(fullfile(pth,[nam '_2mni.mat']),'mni', spm_get_defaults('mat.format'));
    end
    Mat    = mni.affine;
   %do_aff = true;
else
    Mat    = Nii.mat;
   %do_aff = false;
end

% Integrate a Dartel flow field
y0      = spm_dartel_integrate(Nii.dat,job.times,job.K);
if all(job.times == [0 1]),
    mat = Nii.mat0;
    Def = affine(y0,single(Mat));
else
    mat = Mat;
    Def = affine(y0,single(Nii.mat0));
end


%==========================================================================
% function [Def,mat] = get_id(job)
%==========================================================================
function [Def,mat] = get_id(job)
% Get an identity transform based on an image volume
[mat, dim] = spm_get_matdim(job.space{1});
Def = identity(dim, mat);


%==========================================================================
% function [Def,mat] = get_idbbvox(job)
%==========================================================================
function [Def,mat] = get_idbbvox(job)
% Get an identity transform based on bounding box and voxel size.
% This will produce a transversal image.
[mat, dim] = spm_get_matdim('', job.vox, job.bb);
Def = identity(dim, mat);


%==========================================================================
% function [Def,mat] = get_inv(job)
%==========================================================================
function [Def,mat] = get_inv(job)
% Invert a deformation field (derived from a composition of deformations)

NT          = nifti(job.space{:});
[Def0,mat0] = get_comp(job.comp);
M0          = mat0;
M1          = inv(NT.mat);
M0(4,:)     = [0 0 0 1];
M1(4,:)     = [0 0 0 1];
Def         = spm_diffeo('invdef',Def0,NT.dat.dim(1:3),M1,M0);
mat         = NT.mat;
Def         = spm_extrapolate_def(Def,mat);

%==========================================================================
% function fname = save_def(Def,mat,job)
%==========================================================================
function fname = save_def(Def,mat,job)
% Save a deformation field as an image

ofname = job.ofname;
if isempty(ofname), fname = {}; return; end;

[pth,nam] = fileparts(ofname);
if isfield(job.savedir,'savepwd')
    wd = pwd;
elseif isfield(job.savedir,'saveusr')
    wd = job.savedir.saveusr{1};
else
    wd = pwd;
end

fname = {fullfile(wd,['y_' nam '.nii'])};
dim   = [size(Def,1) size(Def,2) size(Def,3) 1 3];
dtype = 'FLOAT32-LE';
off   = 0;
scale = 1;
inter = 0;
dat   = file_array(fname{1},dim,dtype,off,scale,inter);

N      = nifti;
N.dat  = dat;
N.mat  = mat;
N.mat0 = mat;
N.mat_intent  = 'Aligned';
N.mat0_intent = 'Aligned';
N.intent.code = 'VECTOR';
N.intent.name = 'Mapping';
N.descrip     = 'Deformation field';
create(N);
N.dat(:,:,:,1,1) = Def(:,:,:,1);
N.dat(:,:,:,1,2) = Def(:,:,:,2);
N.dat(:,:,:,1,3) = Def(:,:,:,3);


%==========================================================================
% function fname = jac_def(Def,mat,job)
%==========================================================================
function fname = jac_def(Def,mat,job)
% Save Jacobian determinants of deformation field 

ofname = job.ofname;
if isempty(ofname), fname = {}; return; end;

[pth,nam] = fileparts(ofname);
if isfield(job.savedir,'savepwd')
    wd = pwd;
elseif isfield(job.savedir,'saveusr')
    wd = job.savedir.saveusr{1};
else
    wd = pwd;
end

Dets  = spm_diffeo('def2det',Def)/det(mat(1:3,1:3));
Dets(:,:,[1 end]) = NaN;
Dets(:,[1 end],:) = NaN;
Dets([1 end],:,:) = NaN;

fname = {fullfile(wd,['j_' nam '.nii'])};
dim   = [size(Def,1) size(Def,2) size(Def,3) 1 1];
dtype = 'FLOAT32-LE';
off   = 0;
scale = 1;
inter = 0;
dat   = file_array(fname{1},dim,dtype,off,scale,inter);

N      = nifti;
N.dat  = dat;
N.mat  = mat;
N.mat0 = mat;
N.mat_intent  = 'Aligned';
N.mat0_intent = 'Aligned';
N.intent.code = 0;
N.intent.name = 'Mapping';
N.descrip     = 'Jacobian Determinants';
create(N);
N.dat(:,:,:,1,1) = Dets;


%==========================================================================
% function out = pull_def(Def,mat,job)
%==========================================================================
function out = pull_def(Def,mat,job)

PI      = job.fnames;
intrp   = job.interp;
intrp   = [intrp*[1 1 1], 0 0 0];
out     = cell(numel(PI),1);

if numel(PI)==0, return; end

if job.mask
    oM  = zeros(4,4);
    odm = zeros(1,3);
    dim = size(Def);
    msk = true(dim);
    for m=1:numel(PI)
        [pth,nam,ext,num] = spm_fileparts(PI{m});
        NI = nifti(fullfile(pth,[nam ext]));
        dm = NI.dat.dim(1:3);
        if isempty(num)
            j_range = 1:size(NI.dat,4);
        else
            num     = sscanf(num,',%d');
            j_range = num(1);
        end
        for j=j_range

            M0 = NI.mat;
            if ~isempty(NI.extras) && isstruct(NI.extras) && isfield(NI.extras,'mat')
                M1 = NI.extras.mat;
                if size(M1,3) >= j && sum(sum(M1(:,:,j).^2)) ~=0
                    M0 = M1(:,:,j);
                end
            end
            M   = inv(M0);
            if ~all(M(:)==oM(:)) || ~all(dm==odm)
                tmp = affine(Def,M);
                msk = tmp(:,:,:,1)>=1 & tmp(:,:,:,1)<=size(NI.dat,1) ...
                    & tmp(:,:,:,2)>=1 & tmp(:,:,:,2)<=size(NI.dat,2) ...
                    & tmp(:,:,:,3)>=1 & tmp(:,:,:,3)<=size(NI.dat,3);
            end
            oM  = M;
            odm = dm;
        end
    end
end

oM = zeros(4,4);
spm_progress_bar('Init',numel(PI),'Resampling','volumes completed');
for m=1:numel(PI)

    % Generate headers etc for output images
    %----------------------------------------------------------------------
    [pth,nam,ext,num] = spm_fileparts(PI{m});
    NI = nifti(fullfile(pth,[nam ext]));

    j_range = 1:size(NI.dat,4);
    k_range = 1:size(NI.dat,5);
    l_range = 1:size(NI.dat,6);
    if ~isempty(num)
        num = sscanf(num,',%d');
        if numel(num)>=1, j_range = num(1); end
        if numel(num)>=2, k_range = num(2); end
        if numel(num)>=3, l_range = num(3); end
    end

    NO = NI;
    if isfield(job.savedir,'savepwd')
        wd = pwd;
    elseif isfield(job.savedir,'saveusr')
        wd = job.savedir.saveusr{1};
    elseif isfield(job.savedir,'savesrc')
        wd = pth;
    else
        wd = pwd;
    end

    if sum(job.fwhm.^2)==0
        newprefix  = spm_get_defaults('normalise.write.prefix');
        NO.descrip = sprintf('Warped');
    else
        newprefix  = [spm_get_defaults('smooth.prefix') spm_get_defaults('normalise.write.prefix')];
        NO.descrip = sprintf('Smoothed (%gx%gx%g subopt) warped',job.fwhm);
    end
    if isfield(job,'prefix') && ~isempty(job.prefix)
        NO.dat.fname = fullfile(wd,[job.prefix nam ext]);
    else
        NO.dat.fname = fullfile(wd,[newprefix nam ext]);
    end
    dim            = size(Def);
    dim            = dim(1:3);
    NO.dat.dim     = [dim NI.dat.dim(4:end)];
    NO.dat.offset  = 0; % For situations where input .nii images have an extension.
    NO.mat         = mat;
    NO.mat0        = mat;
    NO.mat_intent  = 'Aligned';
    NO.mat0_intent = 'Aligned';
    if isempty(num)
        out{m}     = NO.dat.fname;
    else
        out{m}     = [NO.dat.fname, ',', num2str(num(1))]; 
    end
    NO.extras      = [];
    create(NO);

    % Smoothing settings
    vx  = sqrt(sum(mat(1:3,1:3).^2));
    krn = max(job.fwhm./vx,0.25);

    % Loop over volumes within the file
    %----------------------------------------------------------------------
    %fprintf('%s',nam);
    for j=j_range

        M0 = NI.mat;
        if ~isempty(NI.extras) && isstruct(NI.extras) && isfield(NI.extras,'mat')
            M1 = NI.extras.mat;
            if size(M1,3) >= j && sum(sum(M1(:,:,j).^2)) ~=0
                M0 = M1(:,:,j);
            end
        end
        M  = inv(M0);
        if ~all(M(:)==oM(:))
            % Generate new deformation (if needed)
            Y     = affine(Def,M);
        end
        oM = M;
        % Write the warped data for this time point
        %------------------------------------------------------------------
        for k=k_range
            for l=l_range
                C   = spm_diffeo('bsplinc',single(NI.dat(:,:,:,j,k,l)),intrp);
                dat = spm_diffeo('bsplins',C,Y,intrp);
                if job.mask
                    dat(~msk) = NaN;
                end
                if sum(job.fwhm.^2)~=0
                    spm_smooth(dat,dat,krn); % Side effects
                end
                NO.dat(:,:,:,j,k,l) = dat;
                %fprintf('\t%d,%d,%d', j,k,l);
            end
        end
    end
    %fprintf('\n');
    spm_progress_bar('Set',m);
end
spm_progress_bar('Clear');


%==========================================================================
% function out = push_def(Def,mat,job)
%==========================================================================
function out = push_def(Def,mat,job)
% Generate deformation, which is the inverse of the usual one (it is for "pushing"
% rather than the usual "pulling"). This deformation is affine transformed to
% allow for different voxel sizes and bounding boxes, and also to incorporate
% the affine mapping between MNI space and the population average shape.
%--------------------------------------------------------------------------

% Deal with desired bounding box and voxel sizes.
%--------------------------------------------------------------------------
if isfield(job.fov,'file')
    N1   = nifti(job.fov.file);
    mat0 = N1.mat;
    dim  = N1.dat.dim(1:3);
else
    bb   = job.fov.bbvox.bb;
    vox  = job.fov.bbvox.vox;
    [mat0, dim] = spm_get_matdim('', vox, bb);
end

M   = inv(mat0);
y0  = affine(Def,M);

if isfield(job,'weight') && ~isempty(job.weight) && ~isempty(job.weight{1})
    wfile = job.weight{1};
    Nw    = nifti(wfile);
    Mw    = Nw.mat;
    wt    = Nw.dat(:,:,:,1,1,1);
else
    wt    = [];
end

odm = zeros(1,3);
oM  = zeros(4,4);
PI  = job.fnames;
out = cell(numel(PI),1);
for m=1:numel(PI)

    % Generate headers etc for output images
    %----------------------------------------------------------------------
    [pth,nam,ext,num] = spm_fileparts(PI{m});
    NI = nifti(fullfile(pth,[nam ext]));
    j_range = 1:size(NI.dat,4);
    k_range = 1:size(NI.dat,5);
    l_range = 1:size(NI.dat,6);
    if ~isempty(num)
        num = sscanf(num,',%d');
        if numel(num)>=1, j_range = num(1); end
        if numel(num)>=2, k_range = num(2); end
        if numel(num)>=3, l_range = num(3); end
    end
    
    if isfield(job.savedir,'savepwd')
        wd = pwd;
    elseif isfield(job.savedir,'saveusr')
        wd = job.savedir.saveusr{1};
    elseif isfield(job.savedir,'savesrc')
        wd = pth;
    else
        wd = pwd;
    end

    NO = NI;
    if job.preserve
        NO.dat.scl_slope = 1.0;
        NO.dat.scl_inter = 0.0;
        NO.dat.dtype     = 'float32-le';
        if sum(job.fwhm.^2)==0
            newprefix  = [spm_get_defaults('deformations.modulate.prefix') spm_get_defaults('normalise.write.prefix')];
            NO.descrip = sprintf('Warped & Jac scaled');
        else
            newprefix  = [spm_get_defaults('smooth.prefix') spm_get_defaults('deformations.modulate.prefix') spm_get_defaults('normalise.write.prefix')];
            NO.descrip = sprintf('Smoothed (%gx%gx%g) warped Jac scaled',job.fwhm);
        end
    else
        if sum(job.fwhm.^2)==0
            newprefix  = spm_get_defaults('normalise.write.prefix');
            NO.descrip = sprintf('Warped');
        else
            newprefix  = [spm_get_defaults('smooth.prefix') spm_get_defaults('normalise.write.prefix')];
            NO.descrip = sprintf('Smoothed (%gx%gx%g opt) warped',job.fwhm);
        end
    end
    if isfield(job,'prefix') && ~isempty(job.prefix)
        NO.dat.fname = fullfile(wd,[job.prefix nam ext]);
    else
        NO.dat.fname = fullfile(wd,[newprefix nam ext]);
    end
    NO.dat.dim     = [dim NI.dat.dim(4:end)];
    NO.dat.offset  = 0; % For situations where input .nii images have an extension.
    NO.mat         = mat0;
    NO.mat0        = mat0;
    NO.mat_intent  = 'Aligned';
    NO.mat0_intent = 'Aligned';

    if isempty(num)
        out{m}     = NO.dat.fname;
    else
        out{m}     = [NO.dat.fname, ',', num2str(num(1))];
    end

    NO.extras      = [];
    create(NO);

    % Smoothing settings
    vx  = sqrt(sum(mat0(1:3,1:3).^2));
    krn = max(job.fwhm./vx,0.25);

    % Loop over volumes within the file
    %----------------------------------------------------------------------
    fprintf('%s',nam); drawnow;
    for j=j_range

        % Need to resample the mapping by an affine transform
        % so that it maps from voxels in the native space image
        % to voxels in the spatially normalised image.
        %------------------------------------------------------------------
        M0 = NI.mat;
        if ~isempty(NI.extras) && isstruct(NI.extras) && isfield(NI.extras,'mat')
            M1 = NI.extras.mat;
            if size(M1,3) >= j && sum(sum(M1(:,:,j).^2)) ~=0
                M0 = M1(:,:,j);
            end
        end

        M   = mat\M0;
        dm  = [size(NI.dat),1,1,1,1];
        if ~all(dm(1:3)==odm) || ~all(M(:)==oM(:))
            % Generate new deformation (if needed)
            y   = zeros([dm(1:3),3],'single');
            for d=1:3
                yd = y0(:,:,:,d);
                for x3=1:size(y,3)
                    y(:,:,x3,d) = single(spm_slice_vol(yd,M*spm_matrix([0 0 x3]),dm(1:2),[1 NaN]));
                end
            end
        end

        odm = dm(1:3);
        oM  = M;
        % Write the warped data for this time point.
        %------------------------------------------------------------------
        for k=k_range
            for l=l_range
                f  = single(NI.dat(:,:,:,j,k,l));
                if isempty(wt)
                    if ~job.preserve
                        % Unmodulated - note the slightly novel procedure
                        [f,c] = spm_diffeo('push',f,y,dim);
                        spm_smooth(f,f,krn); % Side effects
                        spm_smooth(c,c,krn); % Side effects
                        f = f./(c+0.001);
                    else
                        % Modulated, by pushing
                        scal = abs(det(NI.mat(1:3,1:3))/det(NO.mat(1:3,1:3))); % Account for vox sizes
                        f    = spm_diffeo('push',f,y,dim)*scal;
                        spm_smooth(f,f,krn); % Side effects
                    end
                else
                    if isequal(size(wt),size(f)) && sum((Mw(:)-M0(:)).^2)<1e-6
                        wtw = single(wt);
                        f   = single(f.*wt);
                    else
                        wtw = zeros(size(f),'single');
                        for z=1:size(wt,3)
                            Mz = Mw\M0*[1 0 0 0; 0 1 0 0; 0 0 1 z; 0 0 0 1];
                            wtw(:,:,z) = single(spm_slice_vol(wt,Mz,[size(f,1),size(f,2)],1));
                        end
                    end
                    if ~job.preserve
                        % Unmodulated - note the slightly novel procedure
                        f = spm_diffeo('push',f.*wtw,y,dim);
                        c = spm_diffeo('push',wtw,y,dim);
                        spm_smooth(f,f,krn); % Side effects
                        spm_smooth(c,c,krn); % Side effects
                        f = f./(c+0.001);
                    else
                        % Modulated, by pushing
                        scal = abs(det(NI.mat(1:3,1:3))/det(NO.mat(1:3,1:3))); % Account for vox sizes
                        f    = spm_diffeo('push',f.*wtw,y,dim)*scal;
                        spm_smooth(f,f,krn); % Side effects
                    end
                    clear wtw
                end
                NO.dat(:,:,:,j,k,l) = f;
                fprintf('\t%d,%d,%d', j,k,l); drawnow;
            end
        end
    end
    fprintf('\n'); drawnow;
end


%==========================================================================
% function out = surf_def(Def,mat,job)
%==========================================================================
function out = surf_def(Def,mat,job)
filenames = job.surface;
out       = cell(numel(filenames),1);
for i=1:numel(filenames)
    fname     = deblank(job.surface{i});
    [pth,nam] = fileparts(fname);
    fprintf('%s\n', nam);
    d      = size(Def);
    tmp    = double(reshape(Def,[d(1:3) 1 d(4)]));
    Tmesh  = spm_swarp(fname, tmp,mat);
    if isfield(job.savedir,'savepwd')
        wd = pwd;
    elseif isfield(job.savedir,'saveusr')
        wd = job.savedir.saveusr{1};
    elseif isfield(job.savedir,'savesrc')
        wd = pth;
    else
        wd = pwd;
    end
    filename = fullfile(wd,[nam,'_warped', '.gii']);
    save(gifti(Tmesh), filename);
    out{i}   = filename;
end


%==========================================================================
% function Def = affine(y,M)
%==========================================================================
function Def = affine(y,M)
Def          = zeros(size(y),'single');
Def(:,:,:,1) = y(:,:,:,1)*M(1,1) + y(:,:,:,2)*M(1,2) + y(:,:,:,3)*M(1,3) + M(1,4);
Def(:,:,:,2) = y(:,:,:,1)*M(2,1) + y(:,:,:,2)*M(2,2) + y(:,:,:,3)*M(2,3) + M(2,4);
Def(:,:,:,3) = y(:,:,:,1)*M(3,1) + y(:,:,:,2)*M(3,2) + y(:,:,:,3)*M(3,3) + M(3,4);


%==========================================================================
% function Def = identity(d,M)
%==========================================================================
function Def = identity(d,M)
[y1,y2]   = ndgrid(single(1:d(1)),single(1:d(2)));
Def       = zeros([d 3],'single');
for y3=1:d(3)
    Def(:,:,y3,1) = y1*M(1,1) + y2*M(1,2) + (y3*M(1,3) + M(1,4));
    Def(:,:,y3,2) = y1*M(2,1) + y2*M(2,2) + (y3*M(2,3) + M(2,4));
    Def(:,:,y3,3) = y1*M(3,1) + y2*M(3,2) + (y3*M(3,3) + M(3,4));
end

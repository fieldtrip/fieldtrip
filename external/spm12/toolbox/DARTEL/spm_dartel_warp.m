function out = spm_dartel_warp(job)
% Register images to template data.
% format spm_dartel_warp(job)
%
% The outputs are flow fields.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dartel_warp.m 4064 2010-09-03 12:57:10Z john $

code = 2;
st = job.settings;
n1 = numel(job.images);
n2 = numel(job.images{1});
NF = struct('NI',[],'vn',[1 1]);
NF(n1,n2) = struct('NI',[],'vn',[1 1]);

for i=1:n1,
    if numel(job.images{i}) ~= n2,
        error('Incompatible number of images');
    end;
    for j=1:n2,
        [pth,nam,ext,num] = spm_fileparts(job.images{i}{j});
        NF(i,j).NI        = nifti(fullfile(pth,[nam ext]));
        num               = [str2num(num) 1 1];
        NF(i,j).vn        = num(1:2);
    end;
end;

dm = [size(NF(1,1).NI.dat) 1];
dm = dm(1:3);

spm_progress_bar('Init',n2,'Registering');
numits = 0;
for it=1:numel(st.param),
    numits = numits + st.param(it).its;
end

for i=1:n2,

    f = zeros([dm,n1],'single');
    for j=1:n1,
        vn = NF(j,i).vn;
        f(:,:,:,j) = single(NF(j,i).NI.dat(:,:,:,vn(1),vn(2)));
    end
    f(~isfinite(f)) = 0;

    [pth,nam,ext] = fileparts(NF(1,i).NI.dat.fname);
    fprintf('*** %s ***\n', nam);

    NU = NF(1,i).NI;
    [pth,nam,ext]   = fileparts(NU.dat.fname);
    NU.dat.fname = fullfile(pth,['u_' nam '.nii']);
    NU.dat.dim   = [dm 1 3];
    NU.dat.dtype = 'float32-le';
    NU.dat.scl_slope = 1;
    NU.dat.scl_inter = 0;
    NU.descrip = 'Flow Field';
    if exist(NU.dat.fname,'file'),
        fprintf('Continuing registration from pre-existing parameters (%s)\n', NU.dat.fname);

        u = NU.dat(:,:,:,1,:);
        u = single(squeeze(u));
    else
        u = zeros([dm,3],'single');
    end;

    it0 = 0;
    for it=1:numel(st.param),
        param = st.param(it);
        prm   = [st.rform, param.rparam, st.optim.lmreg, ...
                 st.optim.cyc, st.optim.its, param.K, code];

        % Load the template
        NG = nifti(strvcat(param.template{:}));

        g  = squeeze(single(NG.dat(:,:,:,:,:)));
        if ~all(size(g)==size(f)),
            error('Incompatible dimensions between images and template');
        end
        for j=1:param.its,
            it0 = it0 + 1;
            [u,ll] = dartel3(u,f,g,prm);
            fprintf('%d \t%g\t%g\t%g\t%g\n',...
                it0,ll(1),ll(2),ll(1)+ll(2),ll(3));

            spm_progress_bar('Set',i-1 + it0/numits);
        end
    end
    create(NU);
    NU.dat(:,:,:,1,:) = reshape(u,[dm 1 3]);
end;
spm_progress_bar('Clear');

n2 = numel(job.images{1});
out.files = cell(n2,1);
for j=1:n2,
    [pth,nam,ext,num] = spm_fileparts(job.images{1}{j});
    fname             = fullfile(pth,['u_' nam '.nii']);
    out.files{j}      = fname;
end;


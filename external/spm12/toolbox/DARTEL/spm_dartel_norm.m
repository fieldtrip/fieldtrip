function out = spm_dartel_norm(job)
% Warp individuals to match template
% FORMAT spm_dartel_norm(job)
% job.flowfields - Flow-fields
% job.images     - Image to warp
% job.interp     - Interpolation method
% job.K          - 2^K timesteps are used
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dartel_norm.m 5668 2013-10-03 18:34:18Z guillaume $


PU = job.flowfields;
PI = job.images;
jactransf = job.jactransf;
intrp    = job.interp;
K        = job.K;

for i=1:numel(PU),
    drawnow;
    NU = nifti(PU{i});
    [pth,nam,ext,num] = spm_fileparts(NU.dat.fname);
    fprintf('%s: ',nam);
    if jactransf,
        [y,dt] = spm_dartel_integrate(NU.dat,[1 0], K);
        dt     = max(dt,eps);
        if jactransf~=1,
            dt = dt.^jactransf;
        end
    else
        y      = spm_dartel_integrate(NU.dat,[1 0], K);
    end;
    y1 = double(y(:,:,:,1));
    y2 = double(y(:,:,:,2));
    y3 = double(y(:,:,:,3));

    for m=1:numel(PI),
        [pth1,nam,ext,num] = spm_fileparts(PI{m}{i});
        NI = nifti(fullfile(pth1,[nam ext]));
        NO = NI;
        if jactransf,
            NO.dat.fname=fullfile(pth,['mw' nam ext]);
            NO.dat.scl_slope = 1.0;
            NO.dat.scl_inter = 0.0;
            NO.dat.dtype     = 'float32-le';
        else
            NO.dat.fname=fullfile(pth,['w' nam ext]);
        end;
        NO.dat.dim = [NU.dat.dim(1:3) NI.dat.dim(4:end)];
        NO.mat  = NU.mat;
        NO.mat0 = NU.mat;
        NO.mat_intent  = 'Aligned';
        NO.mat0_intent = 'Aligned';
        NO.descrip = 'Dartel warped';
        NO.extras  = [];
        create(NO);
        fprintf('%s',nam); drawnow;

        for j=1:size(NI.dat,4),
            if sum(sum((NI.mat  - NU.mat ).^2)) < 0.0001 && ...
               sum(sum((NI.mat0 - NU.mat0).^2)) < 0.0001,
                ty1 = y1;
                ty2 = y2;
                ty3 = y3;
            else
                mat = NI.mat;
                if ~isempty(NI.extras) && isstruct(NI.extras) && isfield(NI.extras,'mat'),
                    mat1 = NI.extras.mat;
                    if size(mat1,3) >= j && sum(sum(mat1(:,:,j).^2)) ~=0,
                        mat = mat1(:,:,j);
                    end;
                end;
                M   = mat\NU.mat0;
                ty1 = M(1,1)*y1 + M(1,2)*y2 + M(1,3)*y3 + M(1,4);
                ty2 = M(2,1)*y1 + M(2,2)*y2 + M(2,3)*y3 + M(2,4);
                ty3 = M(3,1)*y1 + M(3,2)*y2 + M(3,3)*y3 + M(3,4);
            end;
            for k=1:size(NI.dat,5),
                for l=1:size(NI.dat,6),
                    f  = NI.dat(:,:,:,j,k,l);
                    spl_param = [intrp,intrp,intrp,0,0,0];
                    if intrp>1,
                        f = spm_bsplinc(f,spl_param);
                    end;
                    f = spm_bsplins(f,ty1,ty2,ty3,spl_param);
                    if jactransf,
                        NO.dat(:,:,:,j,k,l)=f.*dt;
                    else
                        NO.dat(:,:,:,j,k,l)=f;
                    end;
                    fprintf('\t%d,%d,%d', j,k,l); drawnow;
                end;
            end;
        end;
        fprintf('\n'); drawnow;
    end;
end;

PU = job.flowfields;
PI = job.images;
jactransf = job.jactransf;
out.files = cell(numel(PU),numel(PI));
for i=1:numel(PU),
    [pth,nam] = spm_fileparts(PU{i});
    for m=1:numel(PI),
        [pth1,nam,ext,num] = spm_fileparts(PI{m}{i});
        if jactransf,
            fname = fullfile(pth,['mw' nam ext]);
        else
            fname = fullfile(pth,['w' nam ext]);
        end;
        out.files{i,m} = fname;
    end
end


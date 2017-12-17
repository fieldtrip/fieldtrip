function out = spm_dartel_invnorm(job)
% Warp template to match individuals
% FORMAT spm_dartel_invnorm(job)
% job.flowfields - Filenames of flowfields
% job.images     - Filenames of images to warp
% job.interp     - Interpolation method
% job.K          - 2^K timesteps are used
%
% This function may be useful fo warping labels on to images.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dartel_invnorm.m 5668 2013-10-03 18:34:18Z guillaume $

PU    = job.flowfields;
PI    = job.images;
intrp = job.interp;
K     = job.K;

for i=1:numel(PU),
    [pth1,nam1,ext1,num1] = spm_fileparts(PU{i});
    NU = nifti(fullfile(pth1,[nam1,ext1]));
    fprintf('%s: ',nam1);
    y  = spm_dartel_integrate(NU.dat,[0 1], K);
    y1 = double(y(:,:,:,1));
    y2 = double(y(:,:,:,2));
    y3 = double(y(:,:,:,3));
    clear y

    for m=1:numel(PI),
        [pth2,nam2,ext2,num2] = spm_fileparts(PI{m});
        NI = nifti(fullfile(pth2,[nam2 ext2]));

        NO = NI;
        NO.dat.fname = fullfile(pth1,['w' nam2 '_' nam1 ext2]);
        NO.dat.dim = [NU.dat.dim(1:3) NI.dat.dim(4:end)];
        NO.mat  = NU.mat0;
        NO.mat0 = NU.mat0;
        NO.mat_intent  = NU.mat0_intent;
        NO.mat0_intent = NU.mat0_intent;
        NO.descrip = 'Warped';
        create(NO);
        fprintf('%s',nam2); drawnow;

        for j=1:size(NI.dat,4),
            mat = NI.mat;
            if ~isempty(NI.extras) && isstruct(NI.extras) && isfield(NI.extras,'mat'),
                mat1 = NI.extras.mat;
                if size(mat1,3) >= j && sum(sum(mat1(:,:,j).^2)) ~=0,
                    mat = mat1;
                end;
            end;

            M   = mat\NU.mat;
            ty1 = M(1,1)*y1 + M(1,2)*y2 + M(1,3)*y3 + M(1,4);
            ty2 = M(2,1)*y1 + M(2,2)*y2 + M(2,3)*y3 + M(2,4);
            ty3 = M(3,1)*y1 + M(3,2)*y2 + M(3,3)*y3 + M(3,4);
            for k=1:size(NI.dat,5),
                for l=1:size(NI.dat,6),
                    f             = NI.dat(:,:,:,j,k,l);
                    spl_param     = [intrp,intrp,intrp,0,0,0];
                    if intrp>1, f = spm_bsplinc(f,spl_param); end;
                    f             = spm_bsplins(f,ty1,ty2,ty3,spl_param);
                    NO.dat(:,:,:,j,k,l) = f;
                    fprintf('\t%d,%d,%d', j,k,l); drawnow;
                    clear f
                end;
            end;
            clear ty1 ty2 ty3
        end;
        fprintf('\n'); drawnow;
    end;
end;

PU    = job.flowfields;
PI    = job.images;
out.files = cell(numel(PU),numel(PI));
for i=1:numel(PU),
    [pth1,nam1,ext1,num1] = spm_fileparts(PU{i});
    for m=1:numel(PI),
        [pth2,nam2,ext2,num2] = spm_fileparts(PI{m});
        fname = fullfile(pth1,['w' nam2 '_' nam1 ext2]);
        out.files{i,j} = fname;
    end
end


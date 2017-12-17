function out = spm_dartel_jacobian(job)
% Generate Jacobian determinant fields
% FORMAT spm_dartel_jacobian(job)
% job.flowfields - Filenames of flowfields
% job.K          - 2^K timesteps are used
%
% Note that K needs to be reasonably large in order to obtain reasonable
% Jacobian determinant fields.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dartel_jacobian.m 5506 2013-05-14 17:13:43Z john $

PU = job.flowfields;
K  = job.K;

spm_progress_bar('Init',numel(PU),'Creating Jacobian det fields','Number complete');
for i=1:numel(PU),
    NU = nifti(PU{i});
    [pth,nam,ext,num] = spm_fileparts(NU.dat.fname);
    [y,dt] = spm_dartel_integrate(NU.dat,[1 0], K);
    clear y

    NO = NU;
    NO.dat.fname=fullfile(pth,['jac_' nam(3:end) ext]);
    NO.dat.scl_slope = 1.0;
    NO.dat.scl_inter = 0.0;
    NO.dat.dtype     = 'float32-le';
    NO.dat.dim = NU.dat.dim(1:3);
    NO.mat  = NU.mat;
    NO.mat0 = NU.mat;
    NO.mat_intent  = 'Aligned';
    NO.mat0_intent = 'Aligned';
    NO.descrip = 'Dartel Jacobian';
    create(NO);
    NO.dat(:,:,:)=dt;
    spm_progress_bar('Set',i);
end;
spm_progress_bar('Clear');

PU = job.flowfields;
out.files = cell(numel(PU),1);
for i=1:numel(PU),
    [pth,nam,ext] = fileparts(PU{i});
    fname         = fullfile(pth,['jac_' nam(3:end) ext]);
    out.files{i}  = fname;
end;


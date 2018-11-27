function names = spm_deface(job)
% Face strip images
% FORMAT names = spm_deface(job)
% job.images   - cell array of NIfTI file names
%
% names        - cell array of de-faced images
%
% This is a little routine for attempting to strip the face from images,
% so individuals are more difficult to identify from surface renderings.
%__________________________________________________________________________
% Copyright (C) 2013-2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_deface.m 6086 2014-07-03 16:08:44Z guillaume $ 


if ~nargin
    [P, sts] = spm_select(Inf,'image','Select images to strip the face from');
    if ~sts, return; end
elseif isstruct(job)
    P = job.images;
else
    P = job;
end
P = cellstr(P);
names = cell(size(P));
tpm = spm_load_priors8(spm_vol(fullfile(spm('Dir'),'tpm','TPM.nii')));
for i=1:numel(P)
    names{i} = deface(P{i},tpm);
end

function fname = deface(P,tpm)
nul     = [0 -1.1 0.98 100];
M       = spm_maff8(P,4,20,tpm,[],'mni');
Nii     = nifti(P);
d       = [size(Nii.dat) 1];
[i,j,k] = ndgrid(1:d(1),1:d(2),1:d(3));
nul1    = nul*M*Nii.mat;
msk     = nul1(1)*i + nul1(2)*j + nul1(3)*k + nul1(4) < 0;

fname   = spm_file(Nii.dat.fname,'prefix','anon_');
Noo     = Nii;
Noo.dat.fname = fname;
create(Noo);
for k=1:size(Noo.dat,6),
    for j=1:size(Noo.dat,5),
        for i=1:size(Noo.dat,4), 
            F       = Nii.dat(:,:,:,i,j,k);
            F(msk)  = NaN;
            Noo.dat(:,:,:) = F;
        end
    end
end

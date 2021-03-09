%% Init
[~,ftpath] = ft_version;

%% Read sourcemodel
dat = load(fullfile(ftpath, 'template/sourcemodel/standard_sourcemodel3d10mm.mat'));
sourcemodel = dat.sourcemodel;

%% Normalise #1 - FieldTrip
mri = ft_read_mri(fullfile(ftpath, 'template/anatomy/single_subj_T1.nii'),'dataformat','nifti_spm');
mri.coordsys = 'ras';
cfg     = [];
cfg.dim = mri.dim;
mri     = ft_volumereslice(cfg,mri);

cfg = [];
cfg.spmversion = 'spm12';
cfg.spmmethod = 'new';
cfg.nonlinear = 'no';
normalise = ft_volumenormalise(cfg, mri);

% this will likely crash
pos = ft_warp_apply(normalise.params, sourcemodel.pos, 'sn2individual');

%% Normalise #2 - SPM
% Full version of SPM12 is required

tmpfile = spm_file(tempname,'ext','nii');
copyfile(fullfile(ftpath, 'template/anatomy/single_subj_T1.nii'),tmpfile);
job.channel.vols = {tmpfile};
job.channel.biasreg = 1e-3;
job.channel.biasfwhm = 60;
job.channel.write = [0 0];
job.channel.tpm = fullfile(spm('Dir'),'tpm/TPM.nii');
job.channel.ngaus = [2 2 2 3 4 2];
job.channel.native = [0 0];
job.channel.warped = [0 0];
for t = 1:6
    job.tissue(t).tpm = spm_file(job.channel.tpm,'number',t);
    job.tissue(t).ngaus = job.channel.ngaus(t);
    job.tissue(t).native = [0 0];
    job.tissue(t).warped = [0 0];
end
job.warp.affreg = 'mni';
job.warp.reg = [0 1.0000e-03 0.5000 0.0500 0.2000];
job.warp.samp = 1;
job.warp.write = [0 0];
job.warp.bb = NaN(2,3);
job.warp.vox = 1.5;
job.warp.mrf = 1;
job.warp.cleanup = 0;
spm_preproc_run(job);
seg = load(spm_file(tmpfile,'suffix','_seg8','ext','mat'));

pos = ft_warp_apply(seg, sourcemodel.pos, 'sn2individual');

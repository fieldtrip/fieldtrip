function nmt_sourceoriplot(cfg)
% NMT_SOURCEORIPLOT
% plots functional source orientation data on slices or on
% a surface, optionally as an overlay on anatomical MRI data, where
% statistical data can be used to determine the opacity of the mask. Input
% data comes from FT_SOURCEANALYSIS.
%
% Use as
%   ft_sourceoriplot(cfg, data)
%
% No cfg options are implemented yet.
%
% Requires spm_ov_quivernmt.m to be installed as SPM plugin in
% spm_orthviews subdirectory

global st

% check that spm_ov_quivernmt has been installed in the SPM plug-in directory
if(~exist('spm_orthviews/spm_ov_quivernmt.m','file'))
    error(['Please move or link ' which('spm_ov_quivernmt') ' to: ' fileparts(which('spm_ov_reorient.m'))])
end

if(size(st.nmt.ori,3) > 1)
    ori = st.nmt.ori(:,:,st.nmt.cfg.time_idx(1));
else
    ori = st.nmt.ori;
end
ori = ori ./ sqrt(sum(ori.^2,2)); % normalize

fname{1} = fullfile(tempdir,'oriX.nii');
fname{2} = fullfile(tempdir,'oriY.nii');
fname{3} = fullfile(tempdir,'oriZ.nii');
fname{4} = fullfile(tempdir,'orimask.nii');

for ii=1:3
    oriwrapped{ii} = st.vols{1}.blobs{1}.vol;
    oriwrapped{ii}(st.nmt.cfg.inside_idx) = ori(:,ii);
end
oriwrapped{4} = double(isfinite(oriwrapped{1}));

for ii=1:4
    vol.fname = fname{ii};
    vol.dim = size(st.vols{1}.blobs{1}.vol);
    vol.mat = st.vols{1}.blobs{1}.mat;
    vol.dt = [16 0];
    spm_write_vol(vol,oriwrapped{ii});
end

%%
for ii=1:3
    orivol{ii}=spm_vol(fname{ii});
end
orimaskvol=spm_vol(fname{4});
spm_orthviews('quivernmt','delete',1);
spm_orthviews('quivernmt','init',1,orivol,orimaskvol)
spm_orthviews('redraw')
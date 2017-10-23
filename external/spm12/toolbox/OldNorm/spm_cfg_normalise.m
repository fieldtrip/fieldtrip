function normalise = spm_cfg_normalise
% SPM Configuration file for toolbox 'Old Normalise'
%__________________________________________________________________________
% Copyright (C) 2005-2012 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_normalise.m 4904 2012-09-06 15:08:56Z guillaume $

if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','OldNorm')); end

%--------------------------------------------------------------------------
% source Source Image
%--------------------------------------------------------------------------
source         = cfg_files;
source.tag     = 'source';
source.name    = 'Source Image';
source.help    = {'The image that is warped to match the template(s).  The result is a set of warps, which can be applied to this image, or any other image that is in register with it.'};
source.filter  = 'image';
source.ufilter = '.*';
source.num     = [1 1];

%--------------------------------------------------------------------------
% wtsrc Source Weighting Image
%--------------------------------------------------------------------------
wtsrc         = cfg_files;
wtsrc.tag     = 'wtsrc';
wtsrc.name    = 'Source Weighting Image';
wtsrc.val     = {''};
wtsrc.help    = {'Optional weighting images (consisting of pixel values between the range of zero to one) to be used for registering abnormal or lesioned brains.  These images should match the dimensions of the image from which the parameters are estimated, and should contain zeros corresponding to regions of abnormal tissue.'};
wtsrc.filter  = 'image';
wtsrc.ufilter = '.*';
wtsrc.num     = [0 1];

%--------------------------------------------------------------------------
% subj Subject
%--------------------------------------------------------------------------
subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {source wtsrc };
subj.help    = {'Data for this subject.  The same parameters are used within subject.'};

%--------------------------------------------------------------------------
% esubjs Data
%--------------------------------------------------------------------------
esubjs         = cfg_repeat;
esubjs.tag     = 'esubjs';
esubjs.name    = 'Data';
esubjs.help    = {'List of subjects. Images of each subject should be warped differently.'};
esubjs.values  = {subj };
esubjs.num     = [1 Inf];

%--------------------------------------------------------------------------
% template Template Image
%--------------------------------------------------------------------------
template         = cfg_files;
template.tag     = 'template';
template.name    = 'Template Image';
template.help    = {'Specify a template image to match the source image with. The contrast in the template must be similar to that of the source image in order to achieve a good registration.  It is also possible to select more than one template, in which case the registration algorithm will try to find the best linear combination of these images in order to best model the intensities in the source image.'};
template.filter  = 'image';
template.ufilter = '.*';
template.dir     = fileparts(mfilename('fullpath'));
template.num     = [1 Inf];

%--------------------------------------------------------------------------
% weight Template Weighting Image
%--------------------------------------------------------------------------
weight         = cfg_files;
weight.tag     = 'weight';
weight.name    = 'Template Weighting Image';
weight.val     = {''};
weight.help    = {
                  'Applies a weighting mask to the template(s) during the parameter estimation.  With the default brain mask, weights in and around the brain have values of one whereas those clearly outside the brain are zero.  This is an attempt to base the normalisation purely upon the shape of the brain, rather than the shape of the head (since low frequency basis functions can not really cope with variations in skull thickness).'
                  ''
                  'The option is now available for a user specified weighting image. This should have the same dimensions and mat file as the template images, with values in the range of zero to one.'
}';
weight.filter  = 'image';
weight.ufilter = '.*';
weight.num     = [0 1];

%--------------------------------------------------------------------------
% smosrc Source Image Smoothing
%--------------------------------------------------------------------------
smosrc         = cfg_entry;
smosrc.tag     = 'smosrc';
smosrc.name    = 'Source Image Smoothing';
smosrc.help    = {'Smoothing to apply to a copy of the source image. The template and source images should have approximately the same smoothness. Remember that the templates supplied with SPM have been smoothed by 8mm, and that smoothnesses combine by Pythagoras'' rule.'};
smosrc.strtype = 'e';
smosrc.num     = [1 1];
smosrc.def     = @(val)spm_get_defaults('old.normalise.estimate.smosrc', val{:});

%--------------------------------------------------------------------------
% smoref Template Image Smoothing
%--------------------------------------------------------------------------
smoref         = cfg_entry;
smoref.tag     = 'smoref';
smoref.name    = 'Template Image Smoothing';
smoref.help    = {'Smoothing to apply to a copy of the template image. The template and source images should have approximately the same smoothness. Remember that the templates supplied with SPM have been smoothed by 8mm, and that smoothnesses combine by Pythagoras'' rule.'};
smoref.strtype = 'e';
smoref.num     = [1 1];
smoref.def     = @(val)spm_get_defaults('old.normalise.estimate.smoref', val{:});

%--------------------------------------------------------------------------
% regtype Affine Regularisation
%--------------------------------------------------------------------------
regtype         = cfg_menu;
regtype.tag     = 'regtype';
regtype.name    = 'Affine Regularisation';
regtype.help    = {'Affine registration into a standard space can be made more robust by regularisation (penalising excessive stretching or shrinking).  The best solutions can be obtained by knowing the approximate amount of stretching that is needed (e.g. ICBM templates are slightly bigger than typical brains, so greater zooms are likely to be needed). If registering to an image in ICBM/MNI space, then choose the first option.  If registering to a template that is close in size, then select the second option.  If you do not want to regularise, then choose the third.'};
regtype.labels  = {
                  'ICBM space template'
                  'Average sized template'
                  'No regularisation'
}';
regtype.values  = {
                  'mni'
                  'subj'
                  'none'
}';
regtype.def     = @(val)spm_get_defaults('old.normalise.estimate.regtype', val{:});

%--------------------------------------------------------------------------
% cutoff Nonlinear Frequency Cutoff
%--------------------------------------------------------------------------
cutoff         = cfg_entry;
cutoff.tag     = 'cutoff';
cutoff.name    = 'Nonlinear Frequency Cutoff';
cutoff.help    = {'Cutoff of DCT bases.  Only DCT bases of periods longer than the cutoff are used to describe the warps. The number used will depend on the cutoff and the field of view of the template image(s).'};
cutoff.strtype = 'e';
cutoff.num     = [1 1];
cutoff.def     = @(val)spm_get_defaults('old.normalise.estimate.cutoff', val{:});

%--------------------------------------------------------------------------
% nits Nonlinear Iterations
%--------------------------------------------------------------------------
nits         = cfg_entry;
nits.tag     = 'nits';
nits.name    = 'Nonlinear Iterations';
nits.help    = {'Number of iterations of nonlinear warping performed.'};
nits.strtype = 'w';
nits.num     = [1 1];
nits.def     = @(val)spm_get_defaults('old.normalise.estimate.nits', val{:});

%--------------------------------------------------------------------------
% reg Nonlinear Regularisation
%--------------------------------------------------------------------------
reg         = cfg_entry;
reg.tag     = 'reg';
reg.name    = 'Nonlinear Regularisation';
reg.help    = {'The amount of regularisation for the nonlinear part of the spatial normalisation. Pick a value around one.  However, if your normalised images appear distorted, then it may be an idea to increase the amount of regularisation (by an order of magnitude) - or even just use an affine normalisation. The regularisation influences the smoothness of the deformation fields.'};
reg.strtype = 'e';
reg.num     = [1 1];
reg.def     = @(val)spm_get_defaults('old.normalise.estimate.reg', val{:});

%--------------------------------------------------------------------------
% eoptions Estimation Options
%--------------------------------------------------------------------------
eoptions      = cfg_branch;
eoptions.tag  = 'eoptions';
eoptions.name = 'Estimation Options';
eoptions.val  = {template weight smosrc smoref regtype cutoff nits reg };
eoptions.help = {'Various settings for estimating warps.'};

%--------------------------------------------------------------------------
% est Old Normalise: Estimate
%--------------------------------------------------------------------------
est         = cfg_exbranch;
est.tag     = 'est';
est.name    = 'Old Normalise: Estimate';
est.val     = {esubjs eoptions };
est.help    = {'Computes the warp that best registers a source image (or series of source images) to match a template, saving it to a file imagename''_sn.mat''.'};
est.prog    = @spm_run_normalise;
est.vout    = @vout_estimate;

%--------------------------------------------------------------------------
% matname Parameter File
%--------------------------------------------------------------------------
matname         = cfg_files;
matname.tag     = 'matname';
matname.name    = 'Parameter File';
matname.help    = {'Select the ''_sn.mat'' file containing the spatial normalisation parameters for that subject.'};
matname.filter  = 'mat';
matname.ufilter = '.*_sn\.mat$';
matname.num     = [1 1];

%--------------------------------------------------------------------------
% resample Images to Write
%--------------------------------------------------------------------------
resample         = cfg_files;
resample.tag     = 'resample';
resample.name    = 'Images to Write';
resample.help    = {'These are the images for warping according to the estimated parameters. They can be any images that are in register with the "source" image used to generate the parameters.'};
resample.filter  = 'image';
resample.ufilter = '.*';
resample.num     = [1 Inf];

%--------------------------------------------------------------------------
% subj Subject
%--------------------------------------------------------------------------
subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {matname resample };
subj.help    = {'Data for this subject.  The same parameters are used within subject.'};

%--------------------------------------------------------------------------
% wsubjs Data
%--------------------------------------------------------------------------
wsubjs         = cfg_repeat;
wsubjs.tag     = 'wsubjs';
wsubjs.name    = 'Data';
wsubjs.help    = {'List of subjects. Images of each subject should be warped differently.'};
wsubjs.values  = {subj };
wsubjs.num     = [1 Inf];

%--------------------------------------------------------------------------
% preserve Preserve
%--------------------------------------------------------------------------
preserve         = cfg_menu;
preserve.tag     = 'preserve';
preserve.name    = 'Preserve';
preserve.help    = {
                    'Preserve Concentrations: Spatially normalised images are not "modulated". The warped images preserve the intensities of the original images.'
                    ''
                    'Preserve Total: Spatially normalised images are "modulated" in order to preserve the total amount of signal in the images. Areas that are expanded during warping are correspondingly reduced in intensity.'
}';
preserve.labels = {
                   'Preserve Concentrations'
                   'Preserve Amount'
}';
preserve.values = {0 1};
preserve.def     = @(val)spm_get_defaults('old.normalise.write.preserve', val{:});

%--------------------------------------------------------------------------
% bb Bounding box
%--------------------------------------------------------------------------
bb         = cfg_entry;
bb.tag     = 'bb';
bb.name    = 'Bounding box';
bb.help    = {'The bounding box (in mm) of the volume which is to be written (relative to the anterior commissure).'};
bb.strtype = 'e';
bb.num     = [2 3];
bb.def     = @(val)spm_get_defaults('old.normalise.write.bb', val{:});

%--------------------------------------------------------------------------
% vox Voxel sizes
%--------------------------------------------------------------------------
vox         = cfg_entry;
vox.tag     = 'vox';
vox.name    = 'Voxel sizes';
vox.help    = {'The voxel sizes (x, y & z, in mm) of the written normalised images.'};
vox.strtype = 'e';
vox.num     = [1 3];
vox.def     = @(val)spm_get_defaults('old.normalise.write.vox', val{:});

%--------------------------------------------------------------------------
% interp Interpolation
%--------------------------------------------------------------------------
interp         = cfg_menu;
interp.tag     = 'interp';
interp.name    = 'Interpolation';
interp.help    = {
                  ['The method by which the images are sampled when ' ...
                  'being written in a different space. ' ...
                  '(Note that Inf or NaN values are treated as zero, ' ...
                  'rather than as missing data)']
                  '    Nearest Neighbour:'
                  '      - Fastest, but not normally recommended.'
                  '    Trilinear Interpolation:'
                  '      - OK for PET, realigned fMRI, or segmentations'
                  '    B-spline Interpolation:'
                  ['      - Better quality (but slower) interpolation' ...
                  '/* \cite{thevenaz00a}*/, especially with higher ' ...
                  'degree splines. Can produce values outside the ' ...
                  'original range (e.g. small negative values from an ' ...
                  'originally all positive image).']
}';
interp.labels = {
                 'Nearest neighbour'
                 'Trilinear'
                 '2nd Degree B-spline'
                 '3rd Degree B-Spline '
                 '4th Degree B-Spline '
                 '5th Degree B-Spline'
                 '6th Degree B-Spline'
                 '7th Degree B-Spline'
}';
interp.values = {0 1 2 3 4 5 6 7};
interp.def     = @(val)spm_get_defaults('old.normalise.write.interp', val{:});

%--------------------------------------------------------------------------
% wrap Wrapping
%--------------------------------------------------------------------------
wrap         = cfg_menu;
wrap.tag     = 'wrap';
wrap.name    = 'Wrapping';
wrap.help    = {
                'These are typically:'
                '    No wrapping: for PET or images that have already been spatially transformed. '
                '    Wrap in  Y: for (un-resliced) MRI where phase encoding is in the Y direction (voxel space).'
}';
wrap.labels  = {
               'No wrap'
               'Wrap X'
               'Wrap Y'
               'Wrap X & Y'
               'Wrap Z'
               'Wrap X & Z'
               'Wrap Y & Z'
               'Wrap X, Y & Z'
}';
wrap.values  = {[0 0 0] [1 0 0] [0 1 0] [1 1 0] [0 0 1] [1 0 1] [0 1 1]...
               [1 1 1]};
wrap.def     = @(val)spm_get_defaults('old.normalise.write.wrap', val{:});

%--------------------------------------------------------------------------
% prefix Filename Prefix
%--------------------------------------------------------------------------
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the normalised image file(s). Default prefix is ''w''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.def     = @(val)spm_get_defaults('old.normalise.write.prefix', val{:});

%--------------------------------------------------------------------------
% roptions Writing Options
%--------------------------------------------------------------------------
roptions      = cfg_branch;
roptions.tag  = 'roptions';
roptions.name = 'Writing Options';
roptions.val  = {preserve bb vox interp wrap prefix };
roptions.help = {'Various options for writing normalised images.'};

%--------------------------------------------------------------------------
% write Old Normalise: Write
%--------------------------------------------------------------------------
write         = cfg_exbranch;
write.tag     = 'write';
write.name    = 'Old Normalise: Write';
write.val     = {wsubjs roptions };
write.help    = {'Allows previously estimated warps (stored in imagename''_sn.mat'' files) to be applied to series of images.'};
write.prog    = @spm_run_normalise;
write.vout    = @vout_write;

%--------------------------------------------------------------------------
% source Source Image
%--------------------------------------------------------------------------
source         = cfg_files;
source.tag     = 'source';
source.name    = 'Source Image';
source.help    = {'The image that is warped to match the template(s).  The result is a set of warps, which can be applied to this image, or any other image that is in register with it.'};
source.filter  = 'image';
source.ufilter = '.*';
source.num     = [1 1];

%--------------------------------------------------------------------------
% wtsrc Source Weighting Image
%--------------------------------------------------------------------------
wtsrc         = cfg_files;
wtsrc.tag     = 'wtsrc';
wtsrc.name    = 'Source Weighting Image';
wtsrc.val     = {''};
wtsrc.help    = {'Optional weighting images (consisting of pixel values between the range of zero to one) to be used for registering abnormal or lesioned brains.  These images should match the dimensions of the image from which the parameters are estimated, and should contain zeros corresponding to regions of abnormal tissue.'};
wtsrc.filter  = 'image';
wtsrc.ufilter = '.*';
wtsrc.num     = [0 1];

%--------------------------------------------------------------------------
% resample Images to Write
%--------------------------------------------------------------------------
resample         = cfg_files;
resample.tag     = 'resample';
resample.name    = 'Images to Write';
resample.help    = {'These are the images for warping according to the estimated parameters. They can be any images that are in register with the "source" image used to generate the parameters.'};
resample.filter  = 'image';
resample.ufilter = '.*';
resample.num     = [1 Inf];

%--------------------------------------------------------------------------
% subj Subject
%--------------------------------------------------------------------------
subj         = cfg_branch;
subj.tag     = 'subj';
subj.name    = 'Subject';
subj.val     = {source wtsrc resample };
subj.help    = {'Data for this subject.  The same parameters are used within subject.'};

%--------------------------------------------------------------------------
% ewsubjs Data
%--------------------------------------------------------------------------
ewsubjs         = cfg_repeat;
ewsubjs.tag     = 'ewsubjs';
ewsubjs.name    = 'Data';
ewsubjs.help    = {'List of subjects. Images of each subject should be warped differently.'};
ewsubjs.values  = {subj };
ewsubjs.num     = [1 Inf];

%--------------------------------------------------------------------------
% estwrite Old Normalise: Estimate & Write
%--------------------------------------------------------------------------
estwrite      = cfg_exbranch;
estwrite.tag  = 'estwrite';
estwrite.name = 'Old Normalise: Estimate & Write';
estwrite.val  = {ewsubjs eoptions roptions };
estwrite.help = {'Computes the warp that best registers a source image (or series of source images) to match a template, saving it to the file imagename''_sn.mat''. This option also allows the contents of the imagename''_sn.mat'' files to be applied to a series of images.'};
estwrite.prog = @spm_run_normalise;
estwrite.vout = @vout_estwrite;

%--------------------------------------------------------------------------
% oldnorm Old Normalise
%--------------------------------------------------------------------------
normalise         = cfg_choice;
normalise.tag     = 'oldnorm';
normalise.name    = 'Old Normalise';
normalise.help    = {
                     'This very ancient module spatially (stereotactically) normalises MRI, PET or SPECT images into a standard space defined by some ideal model or template image[s].  The template images supplied with SPM conform to the space defined by the ICBM, NIH P-20 project, and approximate that of the the space described in the atlas of Talairach and Tournoux (1988). The transformation can also be applied to any other image that has been coregistered with these scans. A few researchers may wish to continue using this strategy, but (when good quality anatomical MRI scans are available) the DARTEL approach is now generally recommended instead.'
                     ''
                     'Generally, the algorithms work by minimising the sum of squares difference between the image which is to be normalised, and a linear combination of one or more template images.  For the least squares registration to produce an unbiased estimate of the spatial transformation, the image contrast in the templates (or linear combination of templates) should be similar to that of the image from which the spatial normalisation is derived.  The registration simply searches for an optimum solution.  If the starting estimates are not good, then the optimum it finds may not find the global optimum.'
                     ''
                     'The first step of the normalisation is to determine the optimum 12-parameter affine transformation.  Initially, the registration is performed by matching the whole of the head (including the scalp) to the template.  Following this, the registration proceeded by only matching the brains together, by appropriate weighting of the template voxels.  This is a completely automated procedure (that does not require ``scalp editing'') that discounts the confounding effects of skull and scalp differences.   A Bayesian framework is used, such that the registration searches for the solution that maximises the a posteriori probability of it being correct /* \cite{ashburner97b} */.  i.e., it maximises the product of the likelihood function (derived from the residual squared difference) and the prior function (which is based on the probability of obtaining a particular set of zooms and shears).'
                     ''
                     'The affine registration is followed by estimating nonlinear deformations, whereby the deformations are defined by a linear combination of three dimensional discrete cosine transform (DCT) basis functions /* \cite{ashburner99a} */.  The default options result in each of the deformation fields being described by 1176parameters, where these represent the coefficients of the deformations in three orthogonal directions.  The matching involved simultaneously minimising the membrane energies of the deformation fields and the residual squared difference between the images and template(s).'
                     ''
                     'The primarily use is for stereotactic normalisation to facilitate inter-subject averaging and precise characterisation of functional anatomy /* \cite{ashburner97bir} */.  It is not necessary to spatially normalise the data (this is only a pre-requisite  for  inter-subject averaging or reporting in the Talairach space).  If you wish to circumnavigate this step  (e.g. if you have single slice data or do not have an appropriate high resolution MRI scan) simply specify where you think the  anterior commissure  is  with  the  ORIGIN in the header of the first scan (using the ''Display'' facility) and proceed directly  to ''Smoothing''or ''Statistics''.'
                     ''
                     'All normalised images are written to the same subdirectory as the original images, prefixed with a ''w''.  The details of the transformations are displayed in the results window, and the parameters are saved in the "*_sn.mat" file.'
}';
normalise.values   = {est write estwrite};

%==========================================================================
function dep = vout_estimate(job)
for k=1:numel(job.subj)
    dep(k)            = cfg_dep;
    dep(k).sname      = sprintf('Norm Params File (Subj %d)',k);
    dep(k).src_output = substruct('()',{k},'.','params');
    dep(k).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
end

%==========================================================================
function dep = vout_write(job)
for k=1:numel(job.subj)
    dep(k)            = cfg_dep;
    dep(k).sname      = sprintf('Normalised Images (Subj %d)',k);
    dep(k).src_output = substruct('()',{k},'.','files');
    dep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

%==========================================================================
function dep = vout_estwrite(job)
depe = vout_estimate(job);
depw = vout_write(job);
dep = [depe depw];

function spm_defaults
% Sets the defaults which are used by SPM
%
% FORMAT spm_defaults
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
% Individual users can make copies which can be stored in their own
% matlab subdirectories. If ~/matlab is ahead of the SPM directory
% in the MATLABPATH, then the users own personal defaults are used.
%
% Care must be taken when modifying this file.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner, Andrew Holmes
% $Id$

%-Prevent users from making direct calls to this function
%-----------------------------------------------------------------------
try
    if ~isdeployed
        d = dbstack;
        if isempty(intersect({'spm','spm_get_defaults'},{d.name}))
            fprintf(['Direct calls to spm_defauts are deprecated.\n' ...
                'Please use spm(''defaults'',modality) instead.\n']);
        end
    end
end

global defaults

% Command Line Usage default
%=======================================================================
defaults.cmdline  = 0;

% User Interface defaults
%=======================================================================
defaults.ui.print  = struct('opt',{{'-dpsc2','-append'}},'append',true,'ext','.ps');
defaults.ui.colour = [0.73 0.78 0.96];
defaults.renderer  = 'zbuffer';
defaults.ui.fs     = 14;  % unused

% File format specific
%=======================================================================
% Note that defaults.analyze.flip is no longer used.  Specifying the
% left/right handedness of the voxel indices is now done entirely by
% spm_flip_analyze_images.m
defaults.images.format  = 'img';  % used for DICOM, ECAT and MINC import

% DICOM Import defaults
%=======================================================================
defaults.dicom.root     = 'flat'; % Folder hierarchy

% Stats defaults
%=======================================================================
defaults.stats.maxmem    = 2^26;
defaults.stats.maxres    = 64;
defaults.stats.fmri.ufp  = 0.001;  % Upper tail F-probability
defaults.stats.pet.ufp   = 0.05;
defaults.stats.eeg.ufp   = 1;
defaults.stats.topoFDR   = 1;
defaults.stats.RPVsmooth = 0; % Gaussian FWHM in mm for smoothing RPV image

% Mask defaults
%=======================================================================
defaults.mask.thresh    = 0.8;

% fMRI design defaults
%=======================================================================
defaults.stats.fmri.fmri_t  = 16;
defaults.stats.fmri.fmri_t0 = 1;
defaults.stats.fmri.hpf     = 128;
defaults.stats.fmri.cvi     = 'AR(1)';

% Filename prefix defaults
%=======================================================================
defaults.slicetiming.prefix     = 'a';
defaults.realign.write.prefix   = 'r';
defaults.coreg.write.prefix     = 'r';
defaults.unwarp.write.prefix    = 'u';
defaults.normalise.write.prefix = 'w';
defaults.smooth.prefix          = 's';

% Realignment defaults
%=======================================================================
defaults.realign.estimate.quality = 0.9;
defaults.realign.estimate.weight  = {''};
defaults.realign.estimate.interp  = 2;
defaults.realign.estimate.wrap    = [0 0 0];
defaults.realign.estimate.sep     = 4;
defaults.realign.estimate.fwhm    = 5;
defaults.realign.estimate.rtm     = 1;
defaults.realign.write.mask       = 1;
defaults.realign.write.interp     = 4;
defaults.realign.write.wrap       = [0 0 0];
defaults.realign.write.which      = [2 1];

% Unwarp defaults
%=======================================================================
defaults.unwarp.estimate.rtm      = 0;
defaults.unwarp.estimate.fwhm     = 4;
defaults.unwarp.estimate.basfcn   = [12 12];
defaults.unwarp.estimate.regorder = 1;
defaults.unwarp.estimate.regwgt   = 1e5;
defaults.unwarp.estimate.foe      = [4 5];
defaults.unwarp.estimate.soe      = [];
defaults.unwarp.estimate.rem      = 1;
defaults.unwarp.estimate.jm       = 0;
defaults.unwarp.estimate.noi      = 5;
defaults.unwarp.estimate.expround = 'Average';
%
% Unwarp uses defaults.realign.write 
% defaults for writing.
%

% Coregistration defaults
%=======================================================================
defaults.coreg.estimate.cost_fun = 'nmi';
defaults.coreg.estimate.sep      = [4 2];
defaults.coreg.estimate.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
defaults.coreg.estimate.fwhm     = [7 7];
defaults.coreg.write.interp      = 1;
defaults.coreg.write.wrap        = [0 0 0];
defaults.coreg.write.mask        = 0;

% Spatial Normalisation defaults
%=======================================================================
defaults.normalise.estimate.smosrc  = 8;
defaults.normalise.estimate.smoref  = 0;
defaults.normalise.estimate.regtype = 'mni';
defaults.normalise.estimate.weight  = '';
defaults.normalise.estimate.cutoff  = 25;
defaults.normalise.estimate.nits    = 16;
defaults.normalise.estimate.reg     = 1;
defaults.normalise.estimate.wtsrc   = 0;
defaults.normalise.write.preserve   = 0;
defaults.normalise.write.bb         = [[-78 -112 -50];[78 76 85]];
defaults.normalise.write.vox        = [2 2 2];
defaults.normalise.write.interp     = 1;
defaults.normalise.write.wrap       = [0 0 0];

% VBM Preprocessing defaults
%=======================================================================
defaults.preproc.tpm     = cellstr(char(...
    fullfile(spm('Dir'),'tpm','grey.nii'),...
    fullfile(spm('Dir'),'tpm','white.nii'),...
    fullfile(spm('Dir'),'tpm','csf.nii')));   % Prior probability maps
defaults.preproc.ngaus          = [2 2 2 4]'; % Gaussians per class
defaults.preproc.warpreg        = 1;          % Warping Regularisation
defaults.preproc.warpco         = 25;         % Warp Frequency Cutoff
defaults.preproc.biasreg        = 0.0001;     % Bias regularisation
defaults.preproc.biasfwhm       = 60;         % Bias FWHM
defaults.preproc.regtype        = 'mni';      % Affine Regularisation
defaults.preproc.samp           = 3;          % Sampling distance
defaults.preproc.output.GM      = [0 0 1];
defaults.preproc.output.WM      = [0 0 1];
defaults.preproc.output.CSF     = [0 0 0];
defaults.preproc.output.biascor = 1;
defaults.preproc.output.cleanup = 0;

% Smooth defaults
%=======================================================================
defaults.smooth.fwhm  = [8 8 8];

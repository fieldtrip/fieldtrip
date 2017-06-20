function mri_unbias = ft_volumebiascorrect(cfg, mri)

% FT_VOLUMEBIASCORRECT corrects the image inhomogeneity bias in an anatomical MRI
%
% Use as
%   mri_unbias = ft_volumeunbias(cfg, mri)
%
% See also FT_VOLUMEREALIGN

ft_revision = '$Id$';
ft_nargin = nargin;
ft_nargout = nargout;

% do the general setup of the function

% the ft_preamble function works by calling a number of scripts from
% fieldtrip/utility/private that are able to modify the local workspace

ft_defaults % this ensures that the path is correct and that the ft_defaults global variable is available
ft_preamble init
ft_preamble debug
ft_preamble loadvar     mri
ft_preamble provenance  mri
ft_preamble trackconfig

if ft_abort
  % do not continue function execution in case the outputfile is present and the user indicated to keep it
  return
end

% ensure that the input data is valid for this function
mri = ft_checkdata(mri, 'datatype', {'volume'});

% set the defaults
cfg.spmversion       = ft_getopt(cfg, 'spmversion', 'spm8');
cfg.keepintermediate = ft_getopt(cfg, 'keepintermediate', 'no');

if ~isfield(cfg, 'intermediatename')
  cfg.intermediatename = tempname;
end

% check that the preferred SPM version is on the path
ft_hastoolbox(cfg.spmversion, 1);

% check whether the input has an anatomy
if ~isfield(mri, 'anatomy')
  error('no anatomical information available, this is required for bias correction');
end

% do an approximate alignment
mri_acpc = ft_convert_coordsys(mri, 'spm');

scale = max(mri_acpc.anatomy(:));
mri_acpc.anatomy = mri_acpc.anatomy./scale;

% create an spm-compatible file for the anatomical volume data
V = ft_write_mri([cfg.intermediatename '_anatomy.nii'], mri_acpc.anatomy, 'transform', mri_acpc.transform, 'spmversion', cfg.spmversion, 'dataformat', 'nifti_spm');

switch cfg.spmversion
    case 'spm8'
        
        cfg.nbins = ft_getopt(cfg, 'nbins', 256);
        cfg.reg   = ft_getopt(cfg, 'reg', 0.01);
        cfg.cutoff = ft_getopt(cfg, 'cutoff', 30);
        cfg.niter  = ft_getopt(cfg, 'niter',  128);
        
        flags = struct('nbins', cfg.nbins, 'reg', cfg.reg, 'cutoff', cfg.cutoff, 'niter', cfg.niter);
        T = spm_bias_estimate(V, flags);
        VO = spm_bias_apply(V,T);
        
        mri_unbias = keepfields(mri, {'coordsys', 'dim', 'transform', 'unit', 'inside'});
        mri_unbias.anatomy = VO.dat.*scale;
        
        % if strcmp(cfg.keepintermediate, 'no')
        % % remove the intermediate files
        % for flop=1:length(files)
        % [p, f, x] = fileparts(files{flop});
        % delete(fullfile(p, [f, '.*']));
        % [p, f, x] = fileparts(wfiles{flop});
        % delete(fullfile(p, [f, '.*']));
        % end
        % end
    case 'spm12'
        if ~isfield(cfg, 'tpm') || isempty(cfg.tpm)
            cfg.tpm = fullfile(spm('dir'),'tpm','TPM.nii');
        end
        cfg.biasreg  = ft_getopt(cfg, 'biasreg', 0.0001);
        cfg.biasfwhm = ft_getopt(cfg, 'biasfwhm', 60);
        
        % create the structure that is required for spm_preproc8
        opts.image    = V;
        opts.tpm      = spm_load_priors8(cfg.tpm);
        opts.biasreg  = cfg.biasreg;
        opts.biasfwhm = cfg.biasfwhm;
        opts.lkp      = [1 1 2 2 3 3 4 4 4 5 5 5 5 6 6 ]; % for whatever it's worth
        opts.reg      = [0 0.001 0.5 0.05 0.2];
        opts.samp     = 3;
        opts.fwhm     = 1;
        
        Affine = spm_maff8(opts.image(1),3,32,opts.tpm,eye(4),'mni');
        Affine = spm_maff8(opts.image(1),3, 1,opts.tpm,Affine,'mni');
        opts.Affine = Affine;
        
        % run the segmentation
        fprintf('Estimating the bias field etc..\n');
        p = spm_preproc8(opts);
        
        % this writes the 'native' segmentations
        spm_preproc_write8(p, zeros(6,4), [0 1], [0 1], 0, 0, nan(2,3), nan);
        
        [pathname,name,ext] = fileparts([cfg.intermediatename '_anatomy.nii']);
        filename = fullfile(pathname, ['m' name ext]);
        
        VO = spm_vol(filename);
        dat  = spm_read_vols(VO);
        VO.dat = dat;
        
        mri_unbias = keepfields(mri, {'coordsys', 'dim', 'transform', 'unit', 'inside'});
        mri_unbias.anatomy = VO.dat.*scale;
        
    otherwise
        error('unsupported spmversion requested');
end

ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   mri
ft_postamble provenance mri_unbias
ft_postamble history    mri_unbias
ft_postamble savevar    mri_unbias

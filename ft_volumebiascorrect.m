function mri_unbias = ft_volumebiascorrect(cfg, mri)

% FT_VOLUMEBIASCORRECT corrects the image inhomogeneity bias in an anatomical MRI
%
% Use as
%   mri_unbias = ft_volumebiascorrect(cfg, mri)
% where the input mri should be a single anatomical volume that was for example read with
% FT_READ_MRI. 
%
% The configuration structure can contain
%   cfg.spmversion     = string, 'spm8', 'spm12' (default = 'spm8')
%   cfg.opts           = struct, containing spmversion specific options.
%                        See the code below and the SPM-documentation for
%                        more information.
%
% See also FT_VOLUMEREALIGN FT_VOLUMESEGMENT FT_VOLUMENORMALISE

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
cfg.spmversion       = ft_getopt(cfg, 'spmversion',       'spm12');
cfg.keepintermediate = ft_getopt(cfg, 'keepintermediate', 'no');
cfg.opts             = ft_getopt(cfg, 'opts');

if ~isfield(cfg, 'intermediatename')
  cfg.intermediatename = tempname;
end

% check that the preferred SPM version is on the path
ft_hastoolbox(cfg.spmversion, 1);

% check whether the input has an anatomy
if ~isfield(mri, 'anatomy')
  ft_error('no anatomical information available, this is required for bias correction');
end

% do an approximate alignment
mri_acpc = ft_convert_coordsys(mri, 'acpc');

% flip and permute the 3D volume itself, so that the voxel and
% headcoordinates approximately correspond this improves the convergence
% of the segmentation algorithm
[mri_acpc, permutevec, flipflags] = align_ijk2xyz(mri_acpc);

% create an spm-compatible file for the anatomical volume data
V = ft_write_mri([cfg.intermediatename '_anatomy.nii'], mri_acpc.anatomy, 'transform', mri_acpc.transform, 'spmversion', cfg.spmversion, 'dataformat', 'nifti_spm');

% get the options from hte cfg
opts = cfg.opts;
switch cfg.spmversion
    case 'spm8'
        
        opts.nbins  = ft_getopt(opts, 'nbins', 256);
        opts.reg    = ft_getopt(opts, 'reg', 0.01);
        opts.cutoff = ft_getopt(opts, 'cutoff', 30);
        opts.niter  = ft_getopt(opts, 'niter',  128);
        
        T  = spm_bias_estimate(V, opts);
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
        
        % create the structure that is required for spm_preproc8
        opts.image    = V;
        opts.tpm      = spm_load_priors8(cfg.tpm);
        opts.biasreg  = ft_getopt(opts, 'biasreg',  0.0001);
        opts.biasfwhm = ft_getopt(opts, 'biasfwhm', 60);
        opts.lkp      = ft_getopt(opts, 'lkp',      [1 1 2 2 3 3 4 4 4 5 5 5 5 6 6 ]);
        opts.reg      = ft_getopt(opts, 'reg',      [0 0.001 0.5 0.05 0.2]);
        opts.samp     = ft_getopt(opts, 'samp',     3);
        opts.fwhm     = ft_getopt(opts, 'fwhm',     1);
        
        Affine = spm_maff8(opts.image(1),3,32,opts.tpm,eye(4),'mni');
        Affine = spm_maff8(opts.image(1),3, 1,opts.tpm,Affine,'mni');
        opts.Affine = Affine;
        
        % run the segmentation
        fprintf('Estimating the bias field etc..\n');
        p = spm_preproc8(opts);
        
        % this writes the 'native' segmentations
        spm_preproc_write8(p, [ones(6,1) zeros(6,3)], [0 1], [0 1], 0, 0, nan(2,3), nan);
        
        
        [pathname,name,ext] = fileparts([cfg.intermediatename '_anatomy.nii']);
        filename = fullfile(pathname, ['m' name ext]);
        
        VO     = spm_vol(filename);
        VO.dat = spm_read_vols(VO);
        
        tmp = zeros([size(VO.dat) 6]);
        for k = 1:6
            cfilename    = fullfile(pathname, ['c' num2str(k) name ext]);
            Vtmp         = spm_vol(cfilename);
            tmp(:,:,:,k) = spm_read_vols(Vtmp);
        end
        tmp_tc = false(size(tmp));
        for k = 1:6
            notk = setdiff(1:6,k);
            for m = notk(:)'
                tmp_tc(:,:,:,k) = tmp(:,:,:,m)<tmp(:,:,:,k) & tmp(:,:,:,k);
            end
        end
        diagn = zeros(6,5);
        for k = 1:6
            tmpvals    = tmp(:,:,:,k);
            tmpvals    = tmpvals(tmp_tc(:,:,:,k));
            diagn(k,:) = [sum(sum(sum(tmp_tc(:,:,:,k)))) mean(VO.dat(tmp_tc(:,:,:,k))) std(VO.dat(tmp_tc(:,:,:,k))) mean(tmpvals) std(tmpvals)];
        end
        
   
        mri_unbias = keepfields(mri, {'coordsys', 'dim', 'transform', 'unit', 'inside'});
        mri_unbias.anatomy = VO.dat;
        mri_unbias.params = p;
        mri_unbias.diagn  = diagn;     
    otherwise
        ft_error('unsupported spmversion requested');
end

% flip the volumes back according to the changes introduced by align_ijk2xyz
for k = 1:3
  if flipflags(k)
    mri_unbias.anatomy = flipdim(mri_unbias.anatomy, k);
  end
end

if ~all(permutevec == [1 2 3])
  mri_unbias.anatomy = ipermute(mri_unbias.anatomy, permutevec);
  mri_unbias.dim     = size(mri_unbias.anatomy);
end
    
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   mri
ft_postamble provenance mri_unbias
ft_postamble history    mri_unbias
ft_postamble savevar    mri_unbias

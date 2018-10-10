function [segmented] = ft_volumesegment(cfg, mri)

% FT_VOLUMESEGMENT segments an anatomical MRI. The behaviour depends on the output requested. It can
% return probabilistic tissue maps of gray/white/csf compartments, a skull-stripped anatomy, or
% binary masks representing the brain surface, skull, or scalp surface.
%
% Use as
%   segmented = ft_volumesegment(cfg, mri)
% where the input mri should be a single anatomical volume that was for example read with
% FT_READ_MRI. For the purpose of creating binary masks of the brain or of the skull, you can also
% provide either the anatomical volume or the already segmented volume (with the probabilistic
% tissue maps) as input.
%
% The configuration structure can contain
%   cfg.output         = string or cell-array of strings, see below (default = 'tpm')
%   cfg.spmversion     = string, 'spm2', 'spm8', 'spm12' (default = 'spm12')
%   cfg.spmmethod      = string, 'old', 'new', 'mars' (default = 'old'). This pertains 
%                        to the algorithm used when cfg.spmversion='spm12', see below.
%   cfg.opts           = struct, containing spm-version specific options. See
%                        the code and/or the SPM-documentation for more detail.
%   cfg.template       = filename of the template anatomical MRI (default = '/spm2/templates/T1.mnc'
%                        or '/spm8/templates/T1.nii')
%   cfg.tpm            = cell-array containing the filenames of the tissue probability maps
%   cfg.name           = string for output filename
%   cfg.write          = 'no' or 'yes' (default = 'no'),
%                        writes the probabilistic tissue maps to SPM compatible analyze (spm2),
%                        or nifti (spm8/spm12) files,
%                        with the suffix (spm2)
%                         _seg1, for the gray matter segmentation
%                         _seg2, for the white matter segmentation
%                         _seg3, for the csf segmentation
%                        or with the prefix (spm8, and spm12 with spmmethod='old')
%                         c1, for the gray matter segmentation
%                         c2, for the white matter segmentation
%                         c3, for the csf segmentation
%                        when using spm12 with spmmethod='new' there'll be 3 additional tissue types
%                         c4, for the bone segmentation
%                         c5, for the soft tissue segmentation
%                         c6, for the air segmentation
%                        when using spm12 with spmmethod='mars' the tpms will be
%                         postprocessed with the mars toolbox, yielding smoother%                         segmentations in general.
%   cfg.brainsmooth    = 'no', or scalar, the FWHM of the gaussian kernel in voxels, (default = 5)
%   cfg.scalpsmooth    = 'no', or scalar, the FWHM of the gaussian kernel in voxels, (default = 5)
%   cfg.skullsmooth    = 'no', or scalar, the FWHM of the gaussian kernel in voxels, (default = 5)
%                        this parameter is only used when the segmentation
%                        contains 6 tisuse types, including 'bone'
%   cfg.brainthreshold = 'no', or scalar, relative threshold value which is used to threshold the
%                       tpm in order to create a volumetric brainmask (see below), (default = 0.5)
%   cfg.scalpthreshold = 'no', or scalar, relative threshold value which is used to threshold the
%                        anatomical data in order to create a volumetric scalpmask (see below),
%                        (default = 0.1)
%   cfg.skullthreshold = 'no', or scalar, relative threshold value which is used to threshold the
%                        anatomical data in order to create a volumetric scalpmask (see below),
%                        (default = 0.5). this parameter is only used when
%                        the segmetnation contains 6 tissue types,
%                        including 'bone',
%   cfg.downsample     = integer, amount of downsampling before segmentation
%                       (default = 1; i.e., no downsampling)
%
% The desired segmentation output is specified with cfg.output as a string or cell-array of strings
% and can contain
%   'tpm'         - tissue probability map for csf, white and gray matter
%   'brain'       - binary representation of the brain (including csf, white and gray matter)
%   'skull'       - binary representation of the skull
%   'scalp'       - binary representation of the scalp
%   'skullstrip'  - anatomy with only the brain
%
%
% Example use:
%   cfg        = [];
%   segmented  = ft_volumesegment(cfg, mri) will segmented the anatomy and will output the
%                segmentation result as 3 probabilistic masks in segmented.gray, white and csf.
%
%   cfg.output = 'skullstrip';
%   segmented  = ft_volumesegment(cfg, mri) will generate a skull-stripped anatomy based on a
%                brainmask generated from the probabilistic tissue maps. The skull-stripped anatomy
%                is stored in the field segmented.anatomy.
%
%   cfg.output = {'brain' 'scalp' 'skull'};
%   segmented  = ft_volumesegment(cfg, mri) will produce a volume with 3 binary masks, representing
%                the brain surface, scalp surface, and skull which do not overlap.
%
%   cfg.output = {'scalp'};
%   segmented  = ft_volumesegment(cfg, mri) will produce a volume with a binary mask (based on the
%                anatomy), representing the border of the scalp surface (i.e., everything inside the
%                surface is also included). Such representation of the scalp is produced faster,
%                because it doesn't require to create the tissue probabilty maps before creating
%                the mask.
%
% It is not possible to request tissue-probability map (tpm)  in combination with binary masks
% (brain, scalp or skull) or skull-stripped anatomy. The output will return only the probabilistic
% maps in segmented.gray, white and csf. However, when a segmentation with the probabilistic gray, white
% and csf representations is available, it is possible to use it as input for brain or skull binary mask.
% For example:
%   cfg           = [];
%   cfg.output    = {'tpm'};
%   segment_tpm   = ft_volumesegment(cfg, mri);
%   cfg.output    = {'brain'};
%   segment_brain = ft_volumesegment(cfg, segment_tpm);
%
% For the SPM-based segmentation to work, the coordinate frame of the input MRI needs to be
% approximately coregistered to the templates of the probabilistic tissue maps. The templates are
% defined in SPM/MNI-space. FieldTrip attempts to do an automatic alignment based on the
% coordsys-field in the mri, and if this is not present, based on the coordsys-field in the cfg. If
% none of them is specified the FT_DETERMINE_COORDSYS function is used to interactively assess the
% coordinate system in which the MRI is expressed.
%
% The template mri is defined in SPM/MNI-coordinates:
%   x-axis pointing to the right ear
%   y-axis along the acpc-line
%   z-axis pointing to the top of the head
%   origin in the anterior commissure.
% Note that the segmentation only works if the template MRI is in SPM coordinates.
%
% If the input mri is a string pointing to a CTF *.mri file, the x-axis is assumed to point to the
% nose, and the origin is assumed to be on the interauricular line. In this specific case, when
% ft_read_mri is used to read in the mri, the coordsys field is automatically attached.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat file on disk and/or
% the output data will be written to a *.mat file. These mat files should contain only a single
% variable, corresponding with the input/output structure.
%
% See also FT_READ_MRI, FT_DETERMINE_COORDSYS, FT_PREPARE_HEADMODEL

% Copyright (C) 2007-2017, Jan-Mathijs Schoffelen, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    mri
ft_preamble provenance mri
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% this is not supported any more as of 26/10/2011
if ischar(mri)
  ft_error('please use cfg.inputfile instead of specifying the input variable as a string');
end

% ensure that old and unsupported options are not being relied on by the end-user's script
% instead of specifying cfg.coordsys, the user should specify the coordsys in the data
cfg = ft_checkconfig(cfg, 'forbidden', {'units', 'inputcoordsys', 'coordinates'});
cfg = ft_checkconfig(cfg, 'deprecated',{'coordsys', 'keepintermediate'});
% as of march 2017 keepintermediate is deprecated, does not seem to be
% used, nor sensible. If result files are to be kept, use cfg.write

cfg = ft_checkconfig(cfg, 'renamedval', {'output', 'skin', 'scalp'});

% check if the input data is valid for this function
mri = ft_checkdata(mri, 'datatype', 'volume', 'feedback', 'yes', 'hasunit', 'yes', 'hascoordsys', 'yes');

% set the defaults
cfg.output           = ft_getopt(cfg, 'output',         'tpm');
cfg.downsample       = ft_getopt(cfg, 'downsample',     1);
cfg.spmversion       = ft_getopt(cfg, 'spmversion',     'spm12');
cfg.write            = ft_getopt(cfg, 'write',          'no');
cfg.spmmethod        = ft_getopt(cfg, 'spmmethod',      'old'); % doing old-style in case of spm12

% set default for smooth and threshold
cfg.brainsmooth      = ft_getopt(cfg, 'brainsmooth',    ''); % see also below
cfg.scalpsmooth      = ft_getopt(cfg, 'scalpsmooth',    ''); % see also below
cfg.skullsmooth      = ft_getopt(cfg, 'skullsmooth',    ''); % see also below
cfg.brainthreshold   = ft_getopt(cfg, 'brainthreshold', ''); % see also below
cfg.scalpthreshold   = ft_getopt(cfg, 'scalpthreshold', ''); % see also below

% earlier version of smooth and threshold specification
cfg.smooth           = ft_getopt(cfg, 'smooth',      '');
cfg.threshold        = ft_getopt(cfg, 'threshold',   '');

% chech whether earlier version of smooth and threshold was specified
if ~(isempty(cfg.smooth))
  if isempty(cfg.brainsmooth)
    cfg.brainsmooth = cfg.smooth;
    ft_warning('Smoothing can be specified separately for scalp and brain. User-specified smoothing will be applied for brainmask.')
  end
  if isempty(cfg.scalpsmooth)
    cfg.scalpsmooth = cfg.smooth;
    ft_warning('Smoothing can be specified separately for scalp and brain. User-specified smoothing will be applied for scalpmask.')
  end
end
if ~(isempty(cfg.threshold))
  if isempty(cfg.brainthreshold)
    cfg.brainthreshold = cfg.threshold;
    ft_warning('Threshold can be specified separately for scalp and brain. User-specified threshold will be applied for brainmask.')
  end
  if isempty(cfg.scalpthreshold)
    cfg.scalpthreshold = cfg.threshold;
    ft_warning('Threshold can be specified separately for scalp and brain. User-specified threshold will be applied for scalpmask.')
  end
end
% then set defaults again
cfg.brainsmooth      = ft_getopt(cfg, 'brainsmooth',      5);
cfg.scalpsmooth      = ft_getopt(cfg, 'scalpsmooth',      5);
cfg.skullsmooth      = ft_getopt(cfg, 'skullsmooth',      5); 
cfg.brainthreshold   = ft_getopt(cfg, 'brainthreshold',   0.5);
cfg.scalpthreshold   = ft_getopt(cfg, 'scalpthreshold',   0.1);
cfg.skullthreshold   = ft_getopt(cfg, 'skullthreshold',   0.5);

% check that the preferred SPM version is on the path
ft_hastoolbox(cfg.spmversion, 1);

if ~isfield(cfg, 'name')
  if ~strcmp(cfg.write, 'yes')
    tmp = tempname;
    cfg.name = tmp;
  else
    ft_error('you must specify the output filename in cfg.name');
  end
end
[pathstr, name] = fileparts(cfg.name);
cfg.name = fullfile(pathstr, name); % remove any possible file extension, to be added later

if ~iscell(cfg.output)
  % ensure it to be cell, to allow for multiple outputs
  cfg.output = {cfg.output};
end

% check whether SPM is needed to generate tissue probability maps
if numel(cfg.output) == 1 && strcmp('scalp', cfg.output)
  needtpm = 0; % not needed for (cummulative type) scalpmask
else
  needtpm = any(ismember(cfg.output, {'tpm' 'gray' 'white' 'csf' 'brain' 'skull' 'skullstrip'}));
end

tissue  = isfield(mri, {'gray', 'white', 'csf', 'bone', 'softtissue', 'air'});
ntissue = sum(tissue);
if ntissue==6
  % this is new-style segmentation
  hastpm = ~islogical(mri.gray) && ~islogical(mri.white) && ~islogical(mri.csf);  % tpm should be probabilistic and not binary!
elseif all(tissue(1:3)==true)
  hastpm = ~islogical(mri.gray) && ~islogical(mri.white) && ~islogical(mri.csf);  % tpm should be probabilistic and not binary!
elseif any(tissue)
  ft_warning('the input seems to contain tissue probability maps, but is incomplete for this function');
  hastpm = false;
else
  hastpm = false;
end

if needtpm && ~hastpm
  % spm needs to be used for the creation of the tissue probability maps
  dotpm = 1;
else
  dotpm = 0;
end

if cfg.downsample~=1
  % optionally downsample the anatomical and/or functional volumes
  tmpcfg = keepfields(cfg, {'downsample', 'showcallinfo'});
  tmpcfg.smooth = 'no'; % smoothing is done in ft_volumesegment itself
  mri = ft_volumedownsample(tmpcfg, mri);
  % restore the provenance information
  [cfg, mri] = rollback_provenance(cfg, mri);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the tissue probability maps if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dotpm
  if isdeployed
    % ensure that these exist in this case, otherwise deal with the defaults below
    cfg = ft_checkconfig(cfg, 'required', {'template' 'tpm'});
  end
  
  % remember the original transformation matrix coordinate system
  original = [];
  original.transform = mri.transform;
  original.coordsys  = mri.coordsys;
  if isfield(mri, 'unit')
    original.unit = mri.unit;
  else
    mri = ft_determine_units(mri); % guess the unit field if not present
    original.unit = mri.unit;
  end
  mri = ft_convert_units(mri, 'mm');
  if isdeployed
    mri = ft_convert_coordsys(mri, 'acpc', 2, cfg.template);
  else
    mri = ft_convert_coordsys(mri, 'acpc');
  end

  % flip and permute the 3D volume itself, so that the voxel and
  % headcoordinates approximately correspond this improves the convergence
  % of the segmentation algorithm
  [mri, permutevec, flipflags] = align_ijk2xyz(mri);

  switch lower(cfg.spmversion)
    case 'spm2'
      cfg.template = ft_getopt(cfg, 'template', fullfile(spm('Dir'), 'templates', 'T1.mnc'));
      opts          = ft_getopt(cfg, 'opts');
      opts.estimate = ft_getopt(opts, 'estimate');
      opts.write    = ft_getopt(opts, 'write');
      opts.estimate.affreg = ft_getopt(opts.estimate, 'affreg');      

      VF = ft_write_mri([cfg.name, '.img'], mri.anatomy, 'transform', mri.transform, 'spmversion', cfg.spmversion, 'dataformat', 'analyze_img');

      % set the spm segmentation defaults (from /opt/spm2/spm_defaults.m script)
      opts.estimate.priors = ft_getopt(opts.estimate, 'priors', ...
        char(fullfile(spm('Dir'), 'apriori', 'gray.mnc'),...
             fullfile(spm('Dir'), 'apriori', 'white.mnc'),...
             fullfile(spm('Dir'), 'apriori', 'csf.mnc')));
      opts.estimate.reg    = ft_getopt(opts.estimate, 'reg',    0.01);
      opts.estimate.cutoff = ft_getopt(opts.estimate, 'cutoff', 30);
      opts.estimate.samp   = ft_getopt(opts.estimate, 'samp',   3);
      opts.estimate.bb     = ft_getopt(opts.estimate, 'bb',     [[-88 88]' [-122 86]' [-60 95]']);
      opts.estimate.affreg.smosrc  = ft_getopt(opts.estimate.affreg, 'smosrc',  8);
      opts.estimate.affreg.regtype = ft_getopt(opts.estimate.affreg, 'regtype', 'mni');
      opts.estimate.affreg.weight  = ft_getopt(opts.estimate.affreg, 'weight', '');
      opts.write.cleanup   = ft_getopt(opts.write, 'cleanup', 1);
      opts.write.wrt_cor   = ft_getopt(opts.write, 'wrt_cor', 1);

      % perform the segmentation
      fprintf('performing the segmentation on the specified volume\n');
      spm_segment(VF, cfg.template, opts);
      
      % generate the list of filenames that contains the segmented volumes
      filenames = {[cfg.name, '_seg1.img'];...
                   [cfg.name, '_seg2.img'];...
                   [cfg.name, '_seg3.img']};
      
      
    case 'spm8'
      cfg.tpm = ft_getopt(cfg, 'tpm');
      cfg.tpm = char(cfg.tpm(:));
      if isempty(cfg.tpm)
        cfg.tpm = char(fullfile(spm('Dir'),'tpm','grey.nii'),...
             fullfile(spm('Dir'),'tpm','white.nii'),...
             fullfile(spm('Dir'),'tpm','csf.nii'));
      end
      px.tpm = cfg.tpm;
      
      VF = ft_write_mri([cfg.name, '.img'], mri.anatomy, 'transform', mri.transform, 'spmversion', cfg.spmversion, 'dataformat', 'nifti_spm');

      fprintf('performing the segmentation on the specified volume\n');
      p         = spm_preproc(VF, px);
      [po, dum] = spm_prep2sn(p);
      
      % this writes a mat file, may be needed for Dartel, not sure yet
      save([cfg.name '_sn.mat'],'-struct','po');

      % These settings were taken from a batch
      opts     = ft_getopt(cfg, 'opts');
      opts.GM  = [0 0 1];
      opts.WM  = [0 0 1];
      opts.CSF = [0 0 1];
      opts.biascor = 1;
      opts.cleanup = 0;
      
      % write the segmented volumes, -> this can probably be done differently
      spm_preproc_write(po, opts);

      % generate the list of filenames that contains the segmented volumes
      [pathstr, name] = fileparts(cfg.name);
      filenames = {fullfile(pathstr,['c1', name, '.img']);...
                   fullfile(pathstr,['c2', name, '.img']);...
                   fullfile(pathstr,['c3', name, '.img'])};

    case 'spm12'
      addpath(fullfile(spm('Dir'),'toolbox/OldSeg'));
      if strcmp(cfg.spmmethod, 'old')
        cfg.tpm = ft_getopt(cfg, 'tpm');
        cfg.tpm = char(cfg.tpm(:));
        if isempty(cfg.tpm)
          cfg.tpm = char(fullfile(spm('Dir'),'toolbox/OldSeg','grey.nii'),...
            fullfile(spm('Dir'),'toolbox/OldSeg','white.nii'),...
            fullfile(spm('Dir'),'toolbox/OldSeg','csf.nii'));
        end
        px.tpm = cfg.tpm;
      
        VF = ft_write_mri([cfg.name, '.nii'], mri.anatomy, 'transform', mri.transform, 'spmversion', cfg.spmversion, 'dataformat', 'nifti_spm');
        
        fprintf('performing the segmentation on the specified volume, using the old-style segmentation\n');
        p         = spm_preproc(VF, px);
        [po, dum] = spm_prep2sn(p);
      
        % this write a mat file, may be needed for Dartel, not sure yet
        save([cfg.name '_sn.mat'],'-struct','po');
        
        % These settings are taken from a batch
        opts     = [];
        opts.GM  = [0 0 1];
        opts.WM  = [0 0 1];
        opts.CSF = [0 0 1];
        opts.biascor = 1;
        opts.cleanup = 0;
        
        % write the segmented volumes, -> this can be done differently
        spm_preproc_write(po, opts);
        
        % generate the list of filenames that contains the segmented volumes
        [pathstr, name] = fileparts(cfg.name);
        filenames = {fullfile(pathstr,['c1', name, '.nii']);...
                     fullfile(pathstr,['c2', name, '.nii']);...
                     fullfile(pathstr,['c3', name, '.nii'])};
        
        
      elseif strcmp(cfg.spmmethod, 'new') || strcmp(cfg.spmmethod, 'mars')
        cfg.tpm = ft_getopt(cfg, 'tpm');
        cfg.tpm = char(cfg.tpm);
        if ~isfield(cfg, 'tpm') || isempty(cfg.tpm)
          cfg.tpm = fullfile(spm('dir'),'tpm','TPM.nii');
        end
        
        VF = ft_write_mri([cfg.name, '.nii'], mri.anatomy, 'transform', mri.transform, 'spmversion', cfg.spmversion, 'dataformat', 'nifti_spm');
        
        fprintf('performing the segmentation on the specified volume, using the new-style segmentation\n');
        
        % create the structure that is required for spm_preproc8
        opts          = ft_getopt(cfg, 'opts');
        opts.image    = VF;
        opts.tpm      = ft_getopt(opts, 'tpm', spm_load_priors8(cfg.tpm));
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
        p = spm_preproc8(opts);
        
        % this writes the 'native' segmentations
        if strcmp(cfg.spmmethod, 'new')
          spm_preproc_write8(p, [ones(6,2) zeros(6,2)], [0 0], [0 1], 1, 1, nan(2,3), nan);
        elseif strcmp(cfg.spmmethod, 'mars')
          ft_hastoolbox('mars', 1);
          if ~isfield(cfg, 'mars'), cfg.mars = []; end
          beta        = ft_getopt(cfg.mars, 'beta', 0.1);
          convergence = ft_getopt(cfg.mars, 'convergence', 0.1);
          tcm{1}      = fullfile(fileparts(which('spm_mars_mrf')), 'rTCM_BW20_S1.mat');
          p = spm_mars_mrf(p, [ones(6,2) zeros(6,2)], [0 0], [0 1], tcm, beta, convergence, 1);
        end
        
        % this writes a mat file, may be needed for Dartel, not sure yet
        save([cfg.name '_seg8.mat'],'-struct','p');
        
       
        [pathstr, name] = fileparts(cfg.name);
        filenames = {fullfile(pathstr,['c1', name, '.nii']);...
                     fullfile(pathstr,['c2', name, '.nii']);...
                     fullfile(pathstr,['c3', name, '.nii']);...
                     fullfile(pathstr,['c4', name, '.nii']);...
                     fullfile(pathstr,['c5', name, '.nii']);...
                     fullfile(pathstr,['c6', name, '.nii'])};
        
      else
        ft_error('cfg.spmmethod should be either ''old'', ''new'' or ''mars''');
      end
            
    otherwise
      ft_error('unsupported SPM version');

  end
  
  for k = 1:numel(filenames)
    Vtmp = spm_vol(filenames{k});
    dat  = spm_read_vols(Vtmp);
    Vtmp.dat = dat;
    V(k)     = Vtmp;
  end
  
  if strcmp(cfg.write, 'no')
    [pathstr, name] = fileparts(cfg.name);
    prefix = {'c1';'c2';'c3';'c3';'c4';'c5';'c6';'m';'y_';''};
    suffix = {'_seg1.hdr';'_seg2.hdr';'_seg3.hdr';'_seg1.img';'_seg2.img';'_seg3.img';'_seg1.mat';'_seg2.mat';'_seg3.mat';'.hdr';'.img';'.nii'};
    for k = 1:numel(prefix)
      for j = 1:numel(suffix)
        warning off; % delete is quite verbose in its warnings
        try, delete(fullfile(pathstr,[prefix{k}, name, suffix{j}])); end
        warning on;
      end
    end
  elseif strcmp(cfg.write, 'yes')
    for k= 1:numel(V)
      % I am not sure whether or why this is needed
      V(k).mat = mri.transform;
      V(k)     = spm_create_vol(V(k));
    end
  end
  
  % collect the results
  segmented.dim       = reshape(size(V(1).dat), 1, []); % enforce a row vector
  segmented.transform = original.transform; % use the original transform
  segmented.coordsys  = original.coordsys;  % use the original coordsys
  segmented.unit      = original.unit;      % use the original units
  segmented.gray      = V(1).dat;
  if length(V)>1, segmented.white      = V(2).dat; end
  if length(V)>2, segmented.csf        = V(3).dat; end
  if length(V)>3, segmented.bone       = V(4).dat; end
  if length(V)>4, segmented.softtissue = V(5).dat; end
  if length(V)>5, segmented.air        = V(6).dat; end
  segmented.anatomy   = mri.anatomy;

  % flip the volumes back according to the changes introduced by align_ijk2xyz
  fn = fieldnames(segmented);
  fn = fn(ismember(fn, {'anatomy', 'gray', 'white', 'csf', 'bone', 'softtissue', 'air'}));
  for k = 1:3
    if flipflags(k)
      for j = 1:numel(fn)
        segmented.(fn{j}) = flipdim(segmented.(fn{j}), k);
      end
    end
  end

  if ~all(permutevec == [1 2 3])
    for j = 1:numel(fn)
      segmented.(fn{j}) = ipermute(segmented.(fn{j}), permutevec);
    end
    segmented.dim  = size(segmented.(fn{1}));
  end

else
  % rename the data
  segmented = mri;
  clear mri
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now the data contains the tissue probability maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the requested output fields

remove = {'anatomy' 'csf' 'gray' 'white' 'bone' 'softtissue' 'air' 'skull' 'skullstrip' 'brain' 'scalp'}; % all possible tissues

% check if smoothing or thresholding is required

dosmooth_brain = ~strcmp(cfg.brainsmooth, 'no');
if dosmooth_brain && ischar(cfg.brainsmooth)
  ft_error('invalid value %s for cfg.brainsmooth', cfg.brainsmooth);
end
dosmooth_scalp = ~strcmp(cfg.scalpsmooth, 'no');
if dosmooth_scalp && ischar(cfg.scalpsmooth)
  ft_error('invalid value %s for cfg.scalpsmooth', cfg.scalpsmooth);
end
dosmooth_skull = ~strcmp(cfg.skullsmooth, 'no');
if dosmooth_skull && ischar(cfg.skullsmooth)
  ft_error('invalid value %s for cfg.skullsmooth', cfg.skullsmooth);
end

dothres_brain = ~strcmp(cfg.brainthreshold, 'no');
if dothres_brain && ischar(cfg.brainthreshold)
  ft_error('invalid value %s for cfg.brainthreshold', cfg.brainthreshold);
end
dothres_scalp = ~strcmp(cfg.scalpthreshold, 'no');
if dothres_scalp && ischar(cfg.scalpthreshold)
  ft_error('invalid value %s for cfg.scalpthreshold', cfg.scalpthreshold);
end
dothres_skull = ~strcmp(cfg.skullthreshold, 'no');
if dothres_skull && ischar(cfg.skullthreshold)
  ft_error('invalid value %s for cfg.skullthreshold', cfg.skullthreshold);
end

outp = cfg.output;

if ~isempty(intersect(outp, 'tpm'))
  % output: probability tissue maps
  remove = intersect(remove, {'anatomy'}); 
elseif  ~isempty(intersect(outp, {'white' 'gray' 'csf' 'brain' 'skull' 'scalp' 'skullstrip'}))

  createoutputs = true;
  while createoutputs
    % create scalpmask - no tpm or brainmask is required to create it
    if any(strcmp('scalp', outp))

      fprintf('creating scalpmask ... ');
      % let the softtissue mask take precedence
      if isfield(segmented, 'softtissue')
        fprintf('using the softtissue tpm for segmentation\n');
        fname   = 'softtissue';
        anatomy = segmented.softtissue;
      elseif isfield(segmented, 'anatomy')
        fprintf('using the anatomy field for segmentation\n');
        fname   = 'anatomy';
        anatomy = segmented.anatomy;
      else
        ft_error('no appropriate volume is present for the creation of a scalpmask');
      end
      if dosmooth_scalp
        anatomy = volumesmooth(anatomy, cfg.scalpsmooth, fname);
      else
        fprintf('no smoothing applied on %s for scalp segmentation\n',fname);
      end
      if dothres_scalp
        anatomy = volumethreshold(anatomy, cfg.scalpthreshold, fname);
      else
        fprintf('no threshold applied on %s for scalp segmentation\n',fname);
      end

      % fill the slices along each dimension (because using a single one is
      % just arbitrary, and behavior depends on how the voxeldata is in the
      % volume.
      a1 = volumefillholes(anatomy, 1);
      a2 = volumefillholes(anatomy, 2);
      a3 = volumefillholes(anatomy, 3);

      % anatomy = volumefillholes(anatomy, 2); % FIXME why along the second dimension?
      % scalpmask = anatomy>0;
      scalpmask = a1 | a2 | a3;
      clear anatomy a1 a2 a3

      % threshold again to remove little parts outside of head
      scalpmask = volumethreshold(scalpmask);

      % output: scalp (cummulative) (if this is the only requested output)
      if numel(outp)==1
        segmented.scalp = scalpmask;
		remove(strcmp(remove,'scalp'))=[];
        break
      end
    end   % end scalp

    % create the brain from the tpm
    fprintf('creating brainmask ... using the summation of gray, white and csf tpms\n');
    brain = segmented.gray + segmented.white + segmented.csf;
    if dosmooth_brain
      brain = volumesmooth(brain,  cfg.brainsmooth, 'brainmask');
    else
      fprintf('no smoothing applied on brainmask\n')
    end
    if dothres_brain
      brain = volumethreshold(brain, cfg.brainthreshold, 'brainmask');
    else
      fprintf('no threshold applied on brainmask\n')
    end

    % output: skullstrip
    if any(strcmp('skullstrip', outp))
      if ~isfield(segmented, 'anatomy'), ft_error('no anatomy field present'); end
      fprintf('creating skullstripped anatomy ...');
      brain_ss = cast(brain, class(segmented.anatomy));
      segmented.anatomy = segmented.anatomy.*brain_ss;
      clear brain_ss
      remove(strcmp(remove,'skullstrip'))=[];
	  remove(strcmp(remove,'anatomy'))=[];
      if numel(outp)==1
        break
      end
    end % if skullstrip

    % make binary mask from brain
    brainmask = brain>0;
    clear brain

    % output: brain
    if any(strcmp(outp, 'brain'))
      segmented.brain = brainmask;
	  remove(strcmp(remove,'brain'))=[];
      if numel(outp)==1
        break
      end

      % output: gray, white, csf
    elseif any(strcmp(outp, 'gray')) || any(strcmp(outp, 'white')) || any(strcmp(outp, 'csf'))
      [dum, tissuetype] = max(cat(4, segmented.csf, segmented.gray, segmented.white), [], 4);
      clear dummy
      if any(strcmp(outp, 'white'))
        segmented.white = (tissuetype == 3) & brainmask;
        remove(strcmp(remove,'white'))=[];
      end
      if any(strcmp(outp, 'gray'))
        segmented.gray = (tissuetype == 2) & brainmask;
        remove(strcmp(remove,'gray'))=[];
      end
      if any(strcmp(outp, 'csf'))
        segmented.csf = (tissuetype == 1) & brainmask;
        remove(strcmp(remove,'csf'))=[];
      end

    end % if brain or gray/while/csf

    if any(strcmp('skull', outp)) || any(strcmp('scalp', outp))
      % create skull from brain mask FIXME check this (e.g. strel_bol)
      fprintf('creating skullmask ... ');
      if ~isfield(segmented, 'bone')
        fprintf('using the brainmask\n');
        braindil  = imdilate(brainmask>0, strel_bol(6));
        skullmask = braindil & ~brainmask;
        clear braindil
      elseif isfield(segmented, 'bone')
        fprintf('using the bone tpm for segmentation\n');
        skull = segmented.bone + segmented.gray + segmented.white + segmented.csf;
        if dosmooth_skull
          skull = volumesmooth(skull, cfg.skullsmooth, 'skull');
        else
          fprintf('no smoothing applied on skull tpm for skull segmentation\n');
        end
        if dothres_skull
          skull = volumethreshold(skull, cfg.skullthreshold, 'skull');
        else
          fprintf('no threshold applied on skull tom for skull segmentation\n')
        end
        
        a1 = volumefillholes(skull, 1);
        a2 = volumefillholes(skull, 2);
        a3 = volumefillholes(skull, 3);

        skullmask = a1 | a2 | a3;
        skullmask = volumethreshold(skullmask);
        skullmask = skullmask & ~brainmask;
        
        clear a1 a2 a3
      end
      if any(strcmp(outp, 'skull'))
        segmented.skull = skullmask;
		remove(strcmp(remove,'skull'))=[];
        if numel(outp)==1
          break
        end
      end
      
      % output: scalp (exclusive type)
      if numel(outp) > 1 && any(strcmp('scalp', outp))
        scalpmask(brainmask>0)=0;
        clear brainmask
        scalpmask(skullmask>0)=0;
        clear skullmask
        segmented.scalp=scalpmask;
		remove(strcmp(remove,'scalp'))=[];
        clear scalpmask
      end
    end

    createoutputs = false; % exit the while loop
  end % while createoutputs

else
  ft_error('unknown output %s requested\n', cfg.output{:});
end

% remove unnecessary fields
segmented = removefields(segmented, remove);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous   mri
ft_postamble provenance segmented
ft_postamble history    segmented
ft_postamble savevar    segmented

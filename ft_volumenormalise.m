function [normalised] = ft_volumenormalise(cfg, mri)

% FT_VOLUMENORMALISE normalises anatomical and functional volume data
% to a template anatomical MRI.
%
% Use as
%   [mri] = ft_volumenormalise(cfg, mri)
% where the input mri should be a single anatomical volume that was for
% example read with FT_READ_MRI.
%
% The configuration options can be
%   cfg.parameter        = cell-array with the functional data to be normalised (default = 'all')
%   cfg.keepinside       = 'yes' or 'no', keep the inside/outside labeling (default = 'yes')
%   cfg.downsample       = integer number (default = 1, i.e. no downsampling)
%   cfg.spmversion       = string, 'spm2', 'spm8', 'spm12' (default = 'spm12')
%   cfg.spmmethod        = 'old', 'new' or 'mars', to switch between the different
%                          spm12 implementations. The methods 'new' or 'mars'
%                          uses SPM tissue probability maps instead of the
%                          template MRI specified in cfg.template.
%   cfg.opts             = structure with normalisation options, see SPM documentation for details
%   cfg.template         = string, filename of the template anatomical MRI (default = 'T1.mnc'
%                          for spm2 or 'T1.nii' for spm8 and for spm12).
%   cfg.templatecoordsys = the coordinate system of the template when using a template other
%                          than the default
%   cfg.tpm              = string, file name of the SPM tissue probablility map to use in
%                          case spmversion is 'spm12' and spmmethod is 'new' or 'mars'
%   cfg.write            = 'yes' or 'no' (default = 'no'), writes the segmented volumes to SPM2
%                          compatible analyze-file, with the suffix
%                          _anatomy for the anatomical MRI volume
%                          _param   for each of the functional volumes
%   cfg.name             = string for output filename
%   cfg.keepintermediate = 'yes' or 'no' (default = 'no')
%   cfg.intermediatename = string, prefix of the the coregistered images and of the original
%                          images in the original headcoordinate system
%   cfg.nonlinear        = 'yes' (default) or 'no', estimates a nonlinear transformation
%                          in addition to the linear affine registration. If a reasonably
%                          accurate normalisation is sufficient, a purely linearly transformed
%                          image allows for 'reverse-normalisation', which might come in handy
%                          when for example a region of interest is defined on the normalised
%                          group-average
%   cfg.spmparams        = you can feed in the parameters from a prior normalisation, for example
%                          to apply the parameters determined from an aantomical MRI to an
%                          interpolated source resontruction
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_READ_MRI, FT_VOLUMEDOWNSAMPLE, FT_SOURCEINTERPOLATE, FT_SOURCEPLOT

% Copyright (C) 2004-2020, Jan-Mathijs Schoffelen
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
ft_preamble loadvar mri
ft_preamble provenance mri
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% this is not supported any more as of 26/10/2011
if ischar(mri)
  ft_error('please use cfg.inputfile instead of specifying the input variable as a sting');
end

% ensure that old and unsupported options are not being relied on by the end-user's script
% instead of specifying cfg.coordsys, the user should specify the coordsys in the data
cfg = ft_checkconfig(cfg, 'forbidden', {'units', 'coordsys', 'inputcoord', 'inputcoordsys', 'coordinates'});

% Check coordsys consistency
if isfield(cfg, 'template') && ~isfield(cfg, 'templatecoordsys')
  ft_error('You must specify cfg.templatecoordsys when using your own template')
elseif ~isfield(cfg, 'template') && isfield(cfg, 'templatecoordsys')
  if ~strcmp(cfg.templatecoordsys, 'spm')
    ft_error('Default template must have templatecoordsys=spm')
  end
end

% set the defaults
cfg.spmversion       = ft_getopt(cfg, 'spmversion',       'spm12');
cfg.spmmethod        = ft_getopt(cfg, 'spmmethod',        'old'); % in case of spm12, use the old-style normalisation by default
cfg.opts             = ft_getopt(cfg, 'opts',             []);    % empty will result in default settings
cfg.parameter        = ft_getopt(cfg, 'parameter',        'all');
cfg.downsample       = ft_getopt(cfg, 'downsample',       1);
cfg.write            = ft_getopt(cfg, 'write',            'no');
cfg.keepinside       = ft_getopt(cfg, 'keepinside',       'yes');
cfg.keepintermediate = ft_getopt(cfg, 'keepintermediate', 'no');
cfg.nonlinear        = ft_getopt(cfg, 'nonlinear',        'yes');
cfg.smooth           = ft_getopt(cfg, 'smooth',           'no');
cfg.templatecoordsys = ft_getopt(cfg, 'templatecoordsys', 'spm'); % the assumption is here that the template is one from SPM

% check if the input data is valid for this function
mri = ft_checkdata(mri, 'datatype', 'volume', 'feedback', 'yes', 'hasunit', 'yes', 'hascoordsys', 'yes');

% Warning when using new spmmethod
if isfield(cfg, 'template') && (strcmp(cfg.spmversion, 'spm12') && (strcmp(cfg.spmmethod, 'new') || strcmp(cfg.spmmethod, 'mars')))
  ft_error('spmmethod "%s" only works with the default SPM template', cfg.spmmethod);
end

if ~strcmp(cfg.spmversion, 'spm12') && (strcmp(cfg.spmmethod, 'new') || strcmp(cfg.spmmethod, 'mars'))
  ft_error('spmmethod "%s" only works with SPM version 12', cfg.spmmethod);
end

% check that the preferred SPM version is on the path
ft_hastoolbox(cfg.spmversion, 1);

% check whether the input has an anatomy
if ~isfield(mri, 'anatomy')
  ft_error('no anatomical information available, this is required for normalisation');
end

if isdeployed
  % in deployed mode, FieldTrip cannot use the template in the release version, because these are not compiled
  cfg = ft_checkconfig(cfg, 'required', 'template');
else
  if ~isfield(cfg, 'template')
    spmpath = spm('dir');
    if strcmpi(cfg.spmversion, 'spm2'),  cfg.template = fullfile(spmpath, 'templates', 'T1.mnc'); end
    if strcmpi(cfg.spmversion, 'spm8'),  cfg.template = fullfile(spmpath, 'templates', 'T1.nii'); end
    if strcmpi(cfg.spmversion, 'spm12'), cfg.template = fullfile(spmpath, 'toolbox',   'OldNorm', 'T1.nii'); end
    if ~strcmp(cfg.templatecoordsys, 'spm')
      ft_error('you should specify cfg.templatecoordsys=''spm'' when using an SPM template');
    end
  end
end

% Ensure that the input MRI has interpretable units and that the input MRI is expressed in
% a coordinate system which is in approximate agreement with the template.
ft_notice('Doing initial alignment...')
mri  = ft_convert_units(mri, 'mm');
orig = mri.transform;

mri = ft_convert_coordsys(mri, cfg.templatecoordsys, 2, cfg.template);

% keep track of an initial transformation matrix that does the approximate co-registration
initial = mri.transform / orig;

if strcmp(cfg.keepinside, 'yes')
  % add inside to the list of parameters
  if ~iscell(cfg.parameter)
    cfg.parameter = {cfg.parameter 'inside'};
  else
    cfg.parameter(end+1) = {'inside'};
  end
end

if ~isfield(cfg, 'intermediatename')
  cfg.intermediatename = tempname;
end

if ~isfield(cfg, 'name') && strcmp(cfg.write, 'yes')
  ft_error('you must specify the output filename in cfg.name');
end

if isempty(cfg.template)
  ft_error('you must specify a template anatomical MRI');
end

% the template anatomy should always be stored in a SPM-compatible file
template_ftype = ft_filetype(cfg.template);
if strcmp(template_ftype, 'analyze_hdr') || strcmp(template_ftype, 'analyze_img') || strcmp(template_ftype, 'minc') || strcmp(template_ftype, 'nifti')
  % based on the filetype assume that the coordinates correspond with MNI/SPM convention
  % this is ok
else
  ft_error('the head coordinate system of the template does not seem to be correspond with the mni/spm convention');
end

% select the parameters that should be normalised
cfg.parameter = parameterselection(cfg.parameter, mri);

% the anatomy should always be normalised as the first volume
sel = strcmp(cfg.parameter, 'anatomy');
if ~any(sel)
  cfg.parameter = [{'anatomy'} cfg.parameter];
else
  [dum, indx] = sort(sel);
  cfg.parameter = cfg.parameter(fliplr(indx));
end

if cfg.downsample~=1 && ~(strcmp(cfg.spmversion, 'spm12')&&strcmp(cfg.spmmethod,'new'))
  % optionally downsample the anatomical and/or functional volumes, this is
  % not needed when using spm12 in combination with spmmethod='new'
  tmpcfg = keepfields(cfg, {'downsample', 'parameter', 'smooth', 'showcallinfo'});
  mri = ft_volumedownsample(tmpcfg, mri);
  % restore the provenance information
  [cfg, mri] = rollback_provenance(cfg, mri);
end

% use nii if possible
if strcmpi(cfg.spmversion, 'spm2')
  ext = '.img';
else
  ext = '.nii';
end

% create an spm-compatible header for the anatomical volume data
writeoptions = {'transform',mri.transform,'spmversion',cfg.spmversion};
switch ext
  case '.img'
    % nothing to be done
  case '.nii'
    writeoptions(end+(1:2)) = {'dataformat', 'nifti_spm'};
end
VF = ft_write_mri([cfg.intermediatename '_anatomy' ext], mri.anatomy, writeoptions{:});

% create an spm-compatible file for each of the functional volumes
for k = 2:length(cfg.parameter)  % skip the anatomy
  tmp   = strrep(cfg.parameter{k}, '.', '_');
  data  = reshape(getsubfield(mri, tmp), mri.dim);
  VF(k) = ft_write_mri([cfg.intermediatename '_' tmp ext], data, writeoptions{:});
end

% read the template anatomical volume
switch template_ftype
  case 'minc'
    VG = spm_vol_minc(cfg.template);
  case {'analyze_img', 'analyze_hdr', 'nifti'}
    VG = spm_vol(cfg.template);
  otherwise
    ft_error('Unknown template');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the normalisation parameters, if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(cfg, 'spmparams')
  ft_notice('Performing the parameter estimation');
  
  if strcmp(cfg.spmmethod, 'old') && strcmp(cfg.nonlinear, 'yes')
    ft_info('Warping the individual anatomy to the template anatomy, using non-linear transformations');
    % compute the parameters by warping the individual anatomy
    %VF        = spm_vol([cfg.intermediatename '_anatomy' ext]);
    params    = spm_normalise(VG, VF(1));

  elseif strcmp(cfg.spmmethod, 'old') && strcmp(cfg.nonlinear, 'no')
    ft_info('Warping the individual anatomy to the template anatomy, using only linear transformations');
    % compute the parameters by warping the individual anatomy
    %VF         = spm_vol([cfg.intermediatename '_anatomy' ext]);
    cfg.opts.nits = ft_getopt(cfg.opts, 'nits', 0); % put number of non-linear iterations to zero
    params    = spm_normalise(VG, VF(1), [], [], [], cfg.opts);
  
  elseif strcmp(cfg.spmmethod, 'new') || strcmp(cfg.spmmethod, 'mars')
    ft_info('Warping the individual anatomy to the template anatomy, using the %s-style segmentation', cfg.spmmethod);

    if ~isfield(cfg, 'tpm') || isempty(cfg.tpm)
      ft_notice('Using default tissue probability maps');
      cfg.tpm = fullfile(spm('dir'),'tpm','TPM.nii');
    end
    
    % create the structure that is required for spm_preproc8
    opts          = ft_getopt(cfg, 'opts');
    opts.image    = VF(1);
    opts.tpm      = ft_getopt(opts, 'tpm',      spm_load_priors8(cfg.tpm));
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
    params = spm_preproc8(opts);
    
    % this writes the 'deformation field'
    ft_info('Writing the deformation field to file');
    if strcmp(cfg.spmmethod, 'new')
      bb        = spm_get_bbox(opts.tpm.V(1));
      spm_preproc_write8(params, zeros(6,4), [0 0], [0 1], 1, 1, bb, cfg.downsample);
    elseif strcmp(cfg.spmmethod, 'mars')
      ft_hastoolbox('mars', 1);
      if ~isfield(cfg, 'mars'), cfg.mars = []; end
      beta        = ft_getopt(cfg.mars, 'beta', 0.1);
      convergence = ft_getopt(cfg.mars, 'convergence', 0.1);
      tcm{1}      = fullfile(fileparts(which('spm_mars_mrf')), 'rTCM_BW20_S1.mat');
      params = spm_mars_mrf(params, zeros(6,4), [0 0], [0 1], tcm, beta, convergence, 1);
    else
      ft_error('unsupported spmmethod "%s"', cfg.spmmethod);
    end
    
    oldparams = false;
  end
  
else
  ft_notice('Skipping the parameter estimation, using the parameters specified in the configuration');
  % use the externally specified parameters
  %VF     = spm_vol([cfg.intermediatename '_anatomy' ext]);
  params    = cfg.spmparams;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply the normalisation parameters to the specified volumes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normalised = [];

ft_notice('creating the normalized volumes');
if isfield(params, 'Tr')
  % this is an old-style representation of the parameters
  
  cfg.opts.interp = ft_getopt(cfg.opts, 'interp', 1); % set to 0 for nearest interpolation
  cfg.opts.bb     = ft_getopt(cfg.opts, 'bb',   inf); % set to inf to use template bounding box
  
  % apply the normalisation parameters to each of the volumes
  flags.vox    = cfg.downsample.*[1 1 1];
  flags.interp = cfg.opts.interp;
  flags.bb     = cfg.opts.bb;
  spm_write_sn(char({VF.fname}), params, flags);  % this creates the 'w' prefixed files
  for k = 1:numel(VF)
    [p, f, x] = fileparts(VF(k).fname);
    Vout(k)   = spm_vol(fullfile(p, ['w' f x]));
  end
  
else
  % this is a new- or a mars-style representation of the parameters, it requires spm12 on the path
  ft_hastoolbox('spm12', 1);

  cfg.opts.interp = ft_getopt(cfg.opts, 'interp', 4); % set to 0 for nearest interpolation
  
  [pth,fname,ext] = fileparts(params.image.fname);
  
  tmp        = [];
  tmp.fnames = {VF(:).fname};
  tmp.savedir.saveusr{1} = pth;
  tmp.interp = cfg.opts.interp;
  tmp.mask   = 0;
  tmp.fwhm   = [0 0 0];
  
  job             = [];
  job.comp{1}.def = {fullfile(pth,['y_',fname,ext])};
  job.out{1}.pull = tmp;
  out = spm_deformations(job);
  Vout = spm_vol(char(out.warped));
end

% read the normalised results from the 'w' prefixed files
for k=1:length(Vout)
  normalised = setsubfield(normalised, cfg.parameter{k}, spm_read_vols(Vout(k)));
end

% determine the affine source->template coordinate transformation
final = VG.mat * inv(params.Affine) * inv(VF(1).mat) * initial;

normalised.transform = Vout(1).mat;
normalised.dim       = size(normalised.anatomy);
normalised.params    = params;  % this holds the normalization parameters
normalised.initial   = initial; % this holds the initial co-registration to approximately align with the template
normalised.coordsys  = cfg.templatecoordsys;

if isfield(normalised, 'inside')
  % convert back to a logical volume
  normalised.inside  = abs(normalised.inside-1)<=10*eps;
end

% flip and permute the dimensions to align the volume with the headcoordinate axes
normalised = align_ijk2xyz(normalised);

if strcmp(cfg.write, 'yes')
  % create an spm-compatible file for each of the normalised volumes
  for k = 1:length(cfg.parameter)  % include the anatomy
    tmp  = strrep(cfg.parameter{k}, '.', '_');
    data = reshape(getsubfield(normalised, tmp), normalised.dim);
    ft_write_mri([cfg.name '_' tmp ext], data, writeoptions);
  end
end

if strcmp(cfg.keepintermediate, 'no')
  % remove the intermediate files
  for k = 1:length(Vout)
    [p, f] = fileparts(VF(k).fname);
    delete(fullfile(p, [f, '.*']));
    [p, f] = fileparts(Vout(k).fname);
    delete(fullfile(p, [f, '.*']));
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig

% remember the normalisation parameters in the configuration
% maintain this order for the time being to prevent them to be removed when doing the trackconfig
cfg.spmparams = params;
cfg.final     = final;

ft_postamble previous   mri
ft_postamble provenance normalised
ft_postamble history    normalised
ft_postamble savevar    normalised

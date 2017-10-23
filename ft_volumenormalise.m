function [normalised] = ft_volumenormalise(cfg, mri)

% FT_VOLUMENORMALISE normalises anatomical and functional volume data
% to a template anatomical MRI.
%
% Use as
%   [mri] = ft_volumenormalise(cfg, mri)
% where the input mri should be a single anatomical volume that was for
% example read with FT_READ_MRI.
%
% Configuration options are
%   cfg.spmversion  = string, 'spm2', 'spm8', 'spm12' (default = 'spm8')
%   cfg.template    = string, filename of the template anatomical MRI (default = 'T1.mnc'
%                     for spm2 or 'T1.nii' for spm8)
%   cfg.parameter   = cell-array with the functional data to be normalised (default = 'all')
%   cfg.downsample  = integer number (default = 1, i.e. no downsampling)
%   cfg.name        = string for output filename
%   cfg.write       = 'no' (default) or 'yes', writes the segmented volumes to SPM2
%                     compatible analyze-file, with the suffix
%                     _anatomy for the anatomical MRI volume
%                     _param   for each of the functional volumes
%   cfg.nonlinear   = 'yes' (default) or 'no', estimates a nonlinear transformation
%                     in addition to the linear affine registration. If a reasonably
%                     accurate normalisation is sufficient, a purely linearly transformed
%                     image allows for 'reverse-normalisation', which might come in handy
%                     when for example a region of interest is defined on the normalised
%                     group-average.
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

% Undocumented local options:
%   cfg.keepintermediate = 'yes' or 'no'
%   cfg.intermediatename = prefix of the the coregistered images and of the
%                          original images in the original headcoordinate system
%   cfg.spmparams        = one can feed in parameters from a prior normalisation

% Copyright (C) 2004-2014, Jan-Mathijs Schoffelen
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
cfg = ft_checkconfig(cfg, 'forbidden', {'units', 'inputcoordsys', 'coordinates'});
cfg = ft_checkconfig(cfg, 'deprecated', 'coordsys');
%if isfield(cfg, 'coordsys') && ~isfield(mri, 'coordsys')
%  % from revision 8680 onward (Oct 2013) it is not recommended to use cfg.coordsys to specify the coordinate system of the data.
%  mri.coordsys = cfg.coordsys;
%end

% check if the input data is valid for this function
mri = ft_checkdata(mri, 'datatype', 'volume', 'feedback', 'yes', 'hasunit', 'yes', 'hascoordsys', 'yes');

% set the defaults
cfg.spmversion       = ft_getopt(cfg, 'spmversion',       'spm8');
cfg.spmmethod        = ft_getopt(cfg, 'spmmethod',        'old'); % in case of spm12, use the old-style
cfg.parameter        = ft_getopt(cfg, 'parameter',        'all');
cfg.downsample       = ft_getopt(cfg, 'downsample',       1);
cfg.write            = ft_getopt(cfg, 'write',            'no');
cfg.keepinside       = ft_getopt(cfg, 'keepinside',       'yes');
cfg.keepintermediate = ft_getopt(cfg, 'keepintermediate', 'no');
cfg.nonlinear        = ft_getopt(cfg, 'nonlinear',        'yes');
cfg.smooth           = ft_getopt(cfg, 'smooth',           'no');

% check that the preferred SPM version is on the path
ft_hastoolbox(cfg.spmversion, 1);

% check whether the input has an anatomy
if ~isfield(mri, 'anatomy')
  ft_error('no anatomical information available, this is required for normalisation');
end

% ensure that the data has interpretable units and that the coordinate
% system is in approximate ACPC space and keep track of an initial transformation
% matrix that approximately does the co-registration
mri  = ft_convert_units(mri, 'mm');
orig = mri.transform;
if isdeployed
  mri = ft_convert_coordsys(mri, 'acpc', 2, cfg.template);
else
  mri = ft_convert_coordsys(mri, 'acpc');
end
initial = mri.transform / orig;

if isdeployed
  % in deployed mode, FieldTrip cannot use the template in the release version, because these are not compiled
  cfg = ft_checkconfig(cfg, 'required', 'template');
else
  if ~isfield(cfg, 'template')
    spmpath = spm('dir');
    if strcmpi(cfg.spmversion, 'spm2'),  cfg.template = fullfile(spmpath, filesep, 'templates', filesep, 'T1.mnc'); end
    if strcmpi(cfg.spmversion, 'spm8'),  cfg.template = fullfile(spmpath, filesep, 'templates', filesep, 'T1.nii'); end
    if strcmpi(cfg.spmversion, 'spm12'), cfg.template = fullfile(spmpath, filesep, 'toolbox', filesep, 'OldNorm', filesep, 'T1.nii'); end
  end
end

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

if cfg.downsample~=1
  % optionally downsample the anatomical and/or functional volumes
  tmpcfg = keepfields(cfg, {'downsample', 'parameter', 'smooth', 'showcallinfo'});
  mri = ft_volumedownsample(tmpcfg, mri);
  % restore the provenance information
  [cfg, mri] = rollback_provenance(cfg, mri);
end

ws = ft_warning('off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here the normalisation starts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create an spm-compatible header for the anatomical volume data
VF = ft_write_mri([cfg.intermediatename '_anatomy.img'], mri.anatomy, 'transform', mri.transform, 'spmversion', cfg.spmversion);

% create an spm-compatible file for each of the functional volumes
for parlop=2:length(cfg.parameter)  % skip the anatomy
  tmp  = cfg.parameter{parlop};
  data = reshape(getsubfield(mri, tmp), mri.dim);
  tmp(tmp=='.') = '_';
  ft_write_mri([cfg.intermediatename '_' tmp '.img'], data, 'transform', mri.transform, 'spmversion', cfg.spmversion);
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

fprintf('performing the normalisation\n');
% do spatial normalisation according to these steps
% step 1: read header information for template and source image
% step 2: compute transformation parameters
% step 3: write the results to a file with prefix 'w'

if ~isfield(cfg, 'spmparams') && strcmp(cfg.nonlinear, 'yes')
  fprintf('warping the individual anatomy to the template anatomy\n');
  % compute the parameters by warping the individual anatomy
  VF        = spm_vol([cfg.intermediatename '_anatomy.img']);
  params    = spm_normalise(VG,VF);
elseif ~isfield(cfg, 'spmparams') && strcmp(cfg.nonlinear, 'no')
  fprintf('warping the individual anatomy to the template anatomy, using only linear transformations\n');
  % compute the parameters by warping the individual anatomy
  VF         = spm_vol([cfg.intermediatename '_anatomy.img']);
  flags.nits = 0; % put number of non-linear iterations to zero
  params     = spm_normalise(VG,VF,[],[],[],flags);
else
  fprintf('using the parameters specified in the configuration, skipping the parameter estimation\n');
  % use the externally specified parameters
  VF     = spm_vol([cfg.intermediatename '_anatomy.img']);
  params = cfg.spmparams;
end
flags.vox = [cfg.downsample,cfg.downsample,cfg.downsample];

% determine the affine source->template coordinate transformation
final = VG.mat * inv(params.Affine) * inv(VF.mat) * initial;

% apply the normalisation parameters to each of the volumes
files  = cell(1,numel(cfg.parameter));
wfiles = cell(1,numel(cfg.parameter));
for parlop=1:length(cfg.parameter)
  fprintf('creating normalised analyze-file for %s\n', cfg.parameter{parlop});
  tmp = cfg.parameter{parlop};
  tmp(tmp=='.') = '_';
  files{parlop} = sprintf('%s_%s.img', cfg.intermediatename, tmp);
  [p, f, x] = fileparts(files{parlop});
  wfiles{parlop} = fullfile(p, ['w' f x]);
end
spm_write_sn(char(files),params,flags);  % this creates the 'w' prefixed files

% spm_figure('Create', 'Graphics');
% spm_normalise_disp(params,VF);

normalised = [];

% read the normalised results from the 'w' prefixed files
V = spm_vol(char(wfiles));
for vlop=1:length(V)
  normalised = setsubfield(normalised, cfg.parameter{vlop}, spm_read_vols(V(vlop)));
end

normalised.transform = V(1).mat;
normalised.dim       = size(normalised.anatomy);
normalised.params    = params;  % this holds the normalization parameters
normalised.initial   = initial; % this holds the initial co-registration to approximately align with the template
normalised.coordsys  = 'spm';

if isfield(normalised, 'inside')
  % convert back to a logical volume
  normalised.inside  = abs(normalised.inside-1)<=10*eps;
end

% flip and permute the dimensions to align the volume with the headcoordinate axes
normalised = align_ijk2xyz(normalised);

if strcmp(cfg.write, 'yes')
  % create an spm-compatible file for each of the normalised volumes
  for parlop=1:length(cfg.parameter)  % include the anatomy
    tmp  = cfg.parameter{parlop};
    data = reshape(getsubfield(normalised, tmp), normalised.dim);
    tmp(tmp=='.') = '_';
    ft_write_mri([cfg.name '_' tmp '.img'], data, 'transform', normalised.transform, 'spmversion', cfg.spmversion);
  end
end

if strcmp(cfg.keepintermediate, 'no')
  % remove the intermediate files
  for flop=1:length(files)
    [p, f, x] = fileparts(files{flop});
    delete(fullfile(p, [f, '.*']));
    [p, f, x] = fileparts(wfiles{flop});
    delete(fullfile(p, [f, '.*']));
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig

% remember the normalisation parameters in the configuration
% FIXME maintain this order for the time being to prevent them to be removed when
% doing the trackconfig
cfg.spmparams = params;
cfg.final     = final;

% restore the previous warning state
warning(ws);

ft_postamble previous   mri
ft_postamble provenance normalised
ft_postamble history    normalised
ft_postamble savevar    normalised

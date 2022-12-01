function ft_volumewrite(cfg, volume)

% FT_VOLUMEWRITE exports anatomical or functional volume data to a Analyze
% or BrainVoyager file. The data in the resulting file(s) can be
% further analyzed and/or visualized in MRIcro, SPM, BrainVoyager,
% AFNI or similar packages.
%
% Use as
%   ft_volumewrite(cfg, volume)
% where the input volume structure should represent an anatomical MRI
% that was for example obtained from FT_READ_MRI, the source
% reconstruction results from FT_SOURCEANALYSIS, the statistical
% results from FT_SOURCESTATISTICS or an otherwise processed anatomical
% or functional volume.
%
% The configuration structure should contain the following elements
%   cfg.parameter     = string, describing the functional data to be processed,
%                       e.g. 'pow', 'coh', 'nai' or 'anatomy'
%   cfg.filename      = filename without the extension
%
% To determine the file format, the following option can be specified
%   cfg.filetype      = 'analyze_old', 'nifti' (default), 'nifti_img', 'analyze_spm',
%                       'nifti_spm', 'nifti_gz', 'mgz', 'mgh', 'vmp' or 'vmr'
%
% Depending on the filetype, the cfg should also contain
%   cfg.vmpversion    = 1 or 2, version of the vmp format to use (default = 2)
%   cfg.spmversion    = string, version of SPM to be used (default = 'spm12')
%
% The default filetype is 'nifti', which means that a single *.nii file will be
% written using code from the freesurfer toolbox. The 'nifti_img' filetype uses SPM
% for a dual file (*.img/*.hdr) nifti-format file. The 'nifti_spm' filetype uses SPM
% for a single 'nifti' file.
%
% The analyze, analyze_spm, nifti, nifti_img, nifti_spm and mgz filetypes support a
% homogeneous transformation matrix, the other filetypes do not support a homogeneous
% transformation matrix and hence will be written in their native coordinate system.
%
% You can specify the datatype for the nifti, analyze_spm and analyze_old
% formats. If not specified, the class of the input data will be preserved,
% if the file format allows. Although the higher level function may make an
% attempt to typecast the data, only the nifti fileformat preserves the
% datatype. Also, only when filetype = 'nifti', the slope and intercept
% parameters are stored in the file, so that, when reading the data from
% file, the original values are restored (up to the bit resolution).
%   cfg.datatype      = 'uint8', 'int8', 'int16', 'int32', 'single' or 'double'
%
% By default, integer datatypes will be scaled to the maximum value of the
% physical or statistical parameter, floating point datatypes will not be
% scaled. This can be modified, for instance if the data contains only integers with
% indices into a parcellation, by
%   cfg.scaling       = 'yes' or 'no'
%
% Optional configuration items are
%   cfg.downsample    = integer number (default = 1, i.e. no downsampling)
%   cfg.fiducial.nas  = [x y z] position of nasion
%   cfg.fiducial.lpa  = [x y z] position of LPA
%   cfg.fiducial.rpa  = [x y z] position of RPA
%   cfg.markfiducial  = 'yes' or 'no', mark the fiducials
%   cfg.markorigin    = 'yes' or 'no', mark the origin
%   cfg.markcorner    = 'yes' or 'no', mark the first corner of the volume
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_SOURCEANALYSIS, FT_SOURCESTATISTICS, FT_SOURCEINTERPOLATE, FT_WRITE_MRI

% Undocumented local options:
% cfg.parameter

% Copyright (C) 2003-2006, Robert Oostenveld, Markus Siegel
% Copyright (C) 2011-2021, Jan-Mathijs Schoffelen
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
ft_preamble loadvar volume
ft_preamble provenance volume

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
volume = ft_checkdata(volume, 'datatype', 'volume', 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'forbidden', {'coordsys'});  % the coordinate system should be specified in the data
cfg = ft_checkconfig(cfg, 'required',  {'filename', 'parameter'});
cfg = ft_checkconfig(cfg, 'renamed',   {'coordinates', 'coordsys'});
cfg = ft_checkconfig(cfg, 'renamedval', {'datatype', 'bit1', 'logical'});
cfg = ft_checkconfig(cfg, 'renamedval', {'datatype', 'float', 'single'});
cfg = ft_checkconfig(cfg, 'renamedval', {'filetype', 'analyze', 'analyze_old'});

% set the defaults
cfg.filetype     = ft_getopt(cfg, 'filetype',     'nifti');
cfg.downsample   = ft_getopt(cfg, 'downsample',   1);
cfg.markorigin   = ft_getopt(cfg, 'markorigin',   'no');
cfg.markfiducial = ft_getopt(cfg, 'markfiducial', 'no');
cfg.markcorner   = ft_getopt(cfg, 'markcorner',   'no');
cfg.spmversion   = ft_getopt(cfg, 'spmversion',   'spm12');
cfg.datatype     = ft_getopt(cfg, 'datatype');
cfg.scaling      = ft_getopt(cfg, 'scaling');

if any(strcmp(cfg.datatype, {'logical' 'uint8','int8', 'int16', 'int32'})) && isempty(cfg.scaling)
  cfg.scaling = 'yes';
elseif isempty(cfg.scaling)
  cfg.scaling = 'no';
end

% select the parameter that should be written
cfg.parameter = parameterselection(cfg.parameter, volume);

% only a single parameter should be selected
if iscell(cfg.parameter)
  cfg.parameter = cfg.parameter{1};
end

if cfg.downsample~=1
  % optionally downsample the anatomical and/or functional volumes
  tmpcfg = keepfields(cfg, {'downsample', 'parameter', 'showcallinfo', 'trackcallinfo', 'trackusage', 'trackdatainfo', 'trackmeminfo', 'tracktimeinfo'});
  volume = ft_volumedownsample(tmpcfg, volume);
  % restore the provenance information
  [cfg, volume] = rollback_provenance(cfg, volume);
end

% copy the data and convert into double values so that it can be scaled later
transform = volume.transform;
data      = getsubfield(volume, cfg.parameter);

if strcmp(cfg.markfiducial, 'yes')
  % FIXME THIS DOES NOT WORK, SINCE MINXYZ ETC IS NOT DEFINED
  % FIXME determine the voxel index of the fiducials
  nas = cfg.fiducial.nas;
  lpa = cfg.fiducial.lpa;
  rpa = cfg.fiducial.rpa;
  if any(nas<minxyz) || any(nas>maxxyz)
    ft_warning('nasion does not lie within volume, using nearest voxel');
  end
  if any(lpa<minxyz) || any(lpa>maxxyz)
    ft_warning('LPA does not lie within volume, using nearest voxel');
  end
  if any(rpa<minxyz) || any(rpa>maxxyz)
    ft_warning('RPA does not lie within volume, using nearest voxel');
  end
  idx_nas = [nearest(x, nas(1)) nearest(y, nas(2)) nearest(z, nas(3))];
  idx_lpa = [nearest(x, lpa(1)) nearest(y, lpa(2)) nearest(z, lpa(3))];
  idx_rpa = [nearest(x, rpa(1)) nearest(y, rpa(2)) nearest(z, rpa(3))];
  fprintf('NAS corresponds to voxel [%d, %d, %d]\n', idx_nas);
  fprintf('LPA corresponds to voxel [%d, %d, %d]\n', idx_lpa);
  fprintf('RPA corresponds to voxel [%d, %d, %d]\n', idx_rpa);
  % set the voxel of the fiducials to the maximum value
  data(idx_nas(1), idx_nas(2), idx_nas(3)) = maxval;
  data(idx_lpa(1), idx_lpa(2), idx_lpa(3)) = maxval;
  data(idx_rpa(1), idx_rpa(2), idx_rpa(3)) = maxval;
end

if strcmp(cfg.markorigin, 'yes')
  % FIXME THIS DOES NOT WORK, SINCE MINXYZ ETC IS NOT DEFINED
  % FIXME determine the voxel index of the coordinate system origin
  ori = [0 0 0];
  if any(ori<minxyz) || any(ori>maxxyz)
    ft_warning('origin does not ly within volume, using nearest voxel');
  end
  idx_ori = [nearest(x, ori(1)) nearest(y, ori(2)) nearest(z, ori(3))];
  fprintf('origin corresponds to voxel [%d, %d, %d]\n', idx_ori);
  % set the voxel of the origin to the maximum value
  data(idx_ori(1), idx_ori(2), idx_ori(3)) = maxval;
end

if strcmp(cfg.markcorner, 'yes')
  % set the voxel of the first corner to the maximum value
  data(1:2, 1:1, 1:1) = maxval;     % length 2 along x-axis
  data(1:1, 1:3, 1:1) = maxval;     % length 3 along y-axis
  data(1:1, 1:1, 1:4) = maxval;     % length 4 along z-axis
end

% set not-a-number voxels to zero
data(isnan(data)) = 0;

datatype = class(data);
if isempty(cfg.datatype)
  cfg.datatype = datatype;
end
if ~isequal(datatype, cfg.datatype)
  ft_info('datatype of input data is %s, requested output datatype is %s', datatype, cfg.datatype);
end
if isequal(cfg.datatype, 'logical')
  ft_warning('output datatype of logical is not supported, the data will be stored as uint8');
  cfg.datatype = 'uint8';
end
if istrue(cfg.scaling)
  ft_info('scaling the data and typecasting from %s to %s', datatype, cfg.datatype);
  data   = double(data);
  maxval = max(abs(data(:)));
  minval = min(data(:));
  
  % scale the data so that it fits in the desired numerical data format
  switch lower(cfg.datatype)
    case 'uint8'
      scl_slope = (maxval-minval)/(2^8-1);
      scl_inter = minval;
    case 'int8'
      scl_slope = maxval/(2^7-1);
      scl_inter = 0;
    case 'int16'
      scl_slope = maxval/(2^15-1);
      scl_inter = 0;
    case 'int32'
      scl_slope = maxval/(2^31-1);
      scl_inter = 0;
    case 'single'
      scl_slope = maxval;
      scl_inter = 0;
    case 'double'
      scl_slope = maxval;
      scl_inter = 0;
    otherwise
      ft_error('unknown datatype');
  end
  data = (data - scl_inter)./scl_slope;
else
  % these are parameters that can be written to a nifti file, and can be
  % used to get the data back (close to) its original values
  scl_slope = [];
  scl_inter = [];
end

if ~isequal(datatype, cfg.datatype)
  ft_info('typecasting the numeric data from %s to %s', datatype, cfg.datatype);
  data = cast(data, cfg.datatype);
end

% The BrainVoyager and Analyze format do not support the specification of
% the coordinate system using a homogeneous transformation axis, therefore
% the dimensions of the complete volume has to be reordered by flipping and
% permuting to correspond with their native coordinate system.
switch cfg.filetype
  case {'vmp', 'vmr'}
    if ~isfield(cfg, 'vmpversion')
      fprintf('using BrainVoyager version 2 VMP format\n');
      cfg.vmpversion = 2;
    end

    % the reordering for BrainVoyager has been figured out by Markus Siegel
    if any(strcmp(volume.coordsys, {'als', 'ctf', '4d', 'bti', 'eeglab'}))
      data = permute(data, [2 3 1]);
    elseif any(strcmp(volume.coordsys, {'ras', 'acpc', 'spm', 'mni', 'tal', 'neuromag', 'itab'}))
      data = permute(data, [2 3 1]);
      data = flip(data, 1);
      data = flip(data, 2);
    else
      ft_error('unsupported coordinate system ''%s''', volume.coordsys);
    end
  case 'analyze_old'
    % The coordinate system employed by the ANALYZE programs is left-handed,
    % with the coordinate origin in the lower left corner. Thus, with the
    % subject lying supine, the coordinate origin is on the right side of
    % the body (x), at the back (y), and at the feet (z).
    
    % Analyze   x = right-left
    % Analyze   y = post-ant
    % Analyze   z = inf-sup
    
    % SPM/MNI   x = left-right
    % SPM/MNI   y = post-ant
    % SPM/MNI   z = inf-sup
    
    % CTF       x = post-ant
    % CTF       y = right-left
    % CTF       z = inf-sup

    % the reordering of the Analyze format is according to documentation from Darren Webber
    if any(strcmp(volume.coordsys, {'als', 'ctf', '4d', 'bti', 'eeglab'}))
      data = permute(data, [2 1 3]);
    elseif any(strcmp(volume.coordsys, {'ras', 'acpc', 'spm', 'mni', 'tal', 'neuromag', 'itab'}))
      data = flip(data, 1);
    else
      ft_error('unsupported coordinate system ''%s''', volume.coordsys);
    end
    
  case {'analyze_spm' 'nifti' 'nifti_img' 'nifti_spm' 'nifti_gz' 'mgz' 'mgh'}
    % this format supports a homogenous transformation matrix
    % nothing needs to be changed
  otherwise
    ft_warning('unknown fileformat\n');
end

[p, f, ext] = fileparts(cfg.filename);

% write the volume data to file
switch cfg.filetype
  case {'vmr' 'vmp'}
    ft_write_mri(cfg.filename, data, 'dataformat', cfg.filetype, 'vmpversion', cfg.vmpversion);

  case 'analyze_old'
    if isempty(ext)
      cfg.filename = sprintf('%s.img', cfg.filename);
    end
    ft_write_mri(cfg.filename, data, 'dataformat', 'analyze_old');

  case 'nifti'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write in nifti format, using functions from  the freesurfer toolbox
    % this format supports a homogeneous transformation matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(ext)
      cfg.filename = sprintf('%s.nii', cfg.filename);
    end
    ft_write_mri(cfg.filename, data, 'dataformat', 'nifti', 'transform', transform, 'scl_slope', scl_slope, 'scl_inter', scl_inter);
  
  case 'nifti_gz'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write in zipped nifti format, using functions from  the freesurfer toolbox
    % this format supports a homogeneous transformation matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(ext)
      cfg.filename = sprintf('%s.nii.gz', cfg.filename);
    end
    ft_write_mri(cfg.filename, data, 'dataformat', 'nifti_gz', 'transform', transform, 'scl_slope', scl_slope, 'scl_inter', scl_inter);

  case 'nifti_img'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write in nifti dual file format, using functions from  the SPM toolbox
    % this format supports a homogeneous transformation matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(ext)
      cfg.filename = sprintf('%s.img', cfg.filename);
    end
    ft_write_mri(cfg.filename, data, 'dataformat', 'nifti_img', 'transform', transform, 'spmversion', cfg.spmversion);

  case 'nifti_spm'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write in nifti single file format, using functions from  the SPM toolbox
    % this format supports a homogeneous transformation matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(ext)
      cfg.filename = sprintf('%s.nii', cfg.filename);
    end
    ft_write_mri(cfg.filename, data, 'dataformat', 'nifti_spm', 'transform', transform, 'spmversion', cfg.spmversion);

  case 'analyze_spm'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write in analyze format, using functions from  the SPM toolbox
    % this format supports a homogeneous transformation matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(ext)
      cfg.filename = sprintf('%s.img', cfg.filename);
    end
    ft_write_mri(cfg.filename, data, 'dataformat', 'analyze', 'transform', transform, 'spmversion', cfg.spmversion);

  case {'mgz' 'mgh'}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write in freesurfer_mgz format, using functions from  the freesurfer toolbox
    % this format supports a homogeneous transformation matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ispc && strcmp(cfg.filetype, 'mgz')
      ft_warning('Saving in .mgz format is not possible on a PC, saving in .mgh format instead');
      cfg.filetype = 'mgh';
    end
    if isempty(ext)
      cfg.filename = sprintf('%s.%s', cfg.filename, cfg.filetype);
    end
    ft_write_mri(cfg.filename, data, 'dataformat', cfg.filetype, 'transform', transform);


  otherwise
    fprintf('unknown fileformat\n');
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous volume
ft_postamble provenance

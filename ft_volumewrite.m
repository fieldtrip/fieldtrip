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
%                         e.g. 'pow', 'coh', 'nai' or 'anatomy'
%   cfg.filename      = filename without the extension
%   cfg.filetype      = 'analyze', 'nifti', 'nifti_img', 'analyze_spm', 'mgz',
%                         'vmp' or 'vmr'
%   cfg.vmpversion    = 1 or 2 (default) version of the vmp-format to use
%   cfg.coordsys      = 'spm' or 'ctf', this will only affect the
%                          functionality in case filetype = 'analyze', 'vmp',
%                          or 'vmr'
%
% The default filetype is 'nifti', which means that a single *.nii file
% will be written using the SPM8 toolbox. The 'nifti_img' filetype uses SPM8 for
% a dual file (*.img/*.hdr) nifti-format file.
% The analyze, analyze_spm, nifti, nifti_img and mgz filetypes support a homogeneous
% transformation matrix, the other filetypes do not support a homogeneous transformation
% matrix and hence will be written in their native coordinate system.
%
% You can specify the datatype for the analyze_spm and analyze formats using
%   cfg.datatype      = 'bit1', 'uint8', 'int16', 'int32', 'float' or 'double'
%
% By default, integer datatypes will be scaled to the maximum value of the
% physical or statistical parameter, floating point datatypes will not be
% scaled. This can be modified with
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
% See also FT_SOURCEANALYSIS, FT_SOURCESTATISTICS, FT_SOURCEINTERPOLATE

% Undocumented local options:
% cfg.parameter

% Copyright (C) 2003-2006, Robert Oostenveld, Markus Siegel
% Copyright (C) 2011, Jan-Mathijs Schoffelen
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
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
volume = ft_checkdata(volume, 'datatype', 'volume', 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'required', {'filename', 'parameter'});
cfg = ft_checkconfig(cfg, 'renamed',  {'coordinates', 'coordsys'});

% set the defaults
cfg.filetype     = ft_getopt(cfg, 'filetype',     'nifti');
cfg.datatype     = ft_getopt(cfg, 'datatype',     'int16');
cfg.downsample   = ft_getopt(cfg, 'downsample',   1);
cfg.markorigin   = ft_getopt(cfg, 'markorigin',   'no');
cfg.markfiducial = ft_getopt(cfg, 'markfiducial', 'no');
cfg.markcorner   = ft_getopt(cfg, 'markcorner',   'no');
cfg.scaling      = ft_getopt(cfg, 'scaling',      'no');

if any(strmatch(cfg.datatype, {'int8', 'int16', 'int32'}))
  cfg.scaling = 'yes';
end

if ~isfield(cfg, 'vmpversion') && strcmp(cfg.filetype, 'vmp');
  fprintf('using BrainVoyager version 2 VMP format\n');
  cfg.vmpversion = 2;
end

% select the parameter that should be written
cfg.parameter = parameterselection(cfg.parameter, volume);

% only a single parameter should be selected
if iscell(cfg.parameter)
  cfg.parameter = cfg.parameter{1};
end

if cfg.downsample~=1
  % optionally downsample the anatomical and/or functional volumes
  tmpcfg = keepfields(cfg, {'downsample', 'parameter'});
  volume = ft_volumedownsample(tmpcfg, volume);
  % restore the provenance information
  [cfg, volume] = rollback_provenance(cfg, volume);
end

% copy the data and convert into double values so that it can be scaled later
transform = volume.transform;
data      = double(getsubfield(volume, cfg.parameter));
maxval    = max(data(:));
% ensure that the original volume is not used any more
clear volume

if strcmp(cfg.markfiducial, 'yes')
  % FIXME determine the voxel index of the fiducials
  nas = cfg.fiducial.nas;
  lpa = cfg.fiducial.lpa;
  rpa = cfg.fiducial.rpa;
  if any(nas<minxyz) || any(nas>maxxyz)
    warning('nasion does not lie within volume, using nearest voxel');
  end
  if any(lpa<minxyz) || any(lpa>maxxyz)
    warning('LPA does not lie within volume, using nearest voxel');
  end
  if any(rpa<minxyz) || any(rpa>maxxyz)
    warning('RPA does not lie within volume, using nearest voxel');
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
  % FIXME determine the voxel index of the coordinate system origin
  ori = [0 0 0];
  if any(ori<minxyz) || any(ori>maxxyz)
    warning('origin does not ly within volume, using nearest voxel');
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

if strcmp(cfg.scaling, 'yes')
  % scale the data so that it fits in the desired numerical data format
  switch lower(cfg.datatype)
    case 'bit1'
      data = (data~=0);
    case 'uint8'
      data = uint8((2^8-1) * data./maxval);
    case 'int16'
      data = int16((2^15-1) * data./maxval);
    case 'int32'
      data = int32((2^31-1) * data./maxval);
    case {'single' 'float'}
      data = single(data ./ maxval);
    case 'double'
      data = double(data ./ maxval);
    otherwise
      error('unknown datatype');
  end
end

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

% The BrainVoyager and Analyze format do not support the specification of
% the coordinate system using a homogenous transformation axis, therefore
% the dimensions of the complete volume has to be reordered by flipping and
% permuting to correspond with their native coordinate system.
switch cfg.filetype
  case {'vmp', 'vmr'}
    % the reordering for BrainVoyager has been figured out by Markus Siegel
    if strcmp(cfg.coordsys, 'ctf')
      data = permute(data, [2 3 1]);
    elseif strcmp(cfg.coordsys, 'spm')
      data = permute(data, [2 3 1]);
      data = flipdim(data, 1);
      data = flipdim(data, 2);
    end
    siz = size(data);
  case 'analyze'
    % the reordering of the Analyze format is according to documentation from Darren Webber
    if strcmp(cfg.coordsys, 'ctf')
      data = permute(data, [2 1 3]);
    elseif strcmp(cfg.coordsys, 'spm')
      data = flipdim(data, 1);
    end
    siz = size(data);
  case {'analyze_spm', 'nifti', 'nifti_img' 'mgz' 'mgh'}
    % this format supports a homogenous transformation matrix
    % nothing needs to be changed
  otherwise
    warning('unknown fileformat\n');
end

% write the volume data to file
switch cfg.filetype
  case 'vmp'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write in BrainVoyager VMP format
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen(sprintf('%s.vmp', cfg.filename),'w');
    if fid < 0,
      error('Cannot write to file %s.vmp\n',cfg.filename);
    end

    switch cfg.vmpversion
      case 1
        % write the header
        fwrite(fid, 1, 'short');      % version
        fwrite(fid, 1, 'short');      % number of maps
        fwrite(fid, 1, 'short');      % map type
        fwrite(fid, 0, 'short');      % lag

        fwrite(fid, 0, 'short');      % cluster size
        fwrite(fid, 1, 'float');      % thresh min
        fwrite(fid, maxval, 'float'); % thresh max
        fwrite(fid, 0, 'short');      % df1
        fwrite(fid, 0, 'short');      % df2
        fwrite(fid, 0, 'char');       % name

        fwrite(fid, siz, 'short');    % size
        fwrite(fid, 0, 'short');
        fwrite(fid, siz(1)-1, 'short');
        fwrite(fid, 0, 'short');
        fwrite(fid, siz(2)-1, 'short');
        fwrite(fid, 0, 'short');
        fwrite(fid, siz(3)-1, 'short');
        fwrite(fid, 1, 'short');      % resolution

        % write the data
        fwrite(fid, data, 'float');
      case 2
        % determine relevant subvolume
        % FIXME, this is not functional at the moment, since earlier in this function all nans have been replaced by zeros
        minx = min(find(~isnan(max(max(data,[],3),[],2))));
        maxx = max(find(~isnan(max(max(data,[],3),[],2))));
        miny = min(find(~isnan(max(max(data,[],3),[],1))));
        maxy = max(find(~isnan(max(max(data,[],3),[],1))));
        minz = min(find(~isnan(max(max(data,[],1),[],2))));
        maxz = max(find(~isnan(max(max(data,[],1),[],2))));

        % write the header
        fwrite(fid, 2, 'short');      % version
        fwrite(fid, 1, 'int');        % number of maps
        fwrite(fid, 1, 'int');        % map type
        fwrite(fid, 0, 'int');        % lag

        fwrite(fid, 0, 'int');        % cluster size
        fwrite(fid, 0, 'char');       % cluster enable
        fwrite(fid, 1, 'float');      % thresh
        fwrite(fid, maxval, 'float'); % thresh
        fwrite(fid, 0, 'int');        % df1
        fwrite(fid, 0, 'int');        % df2
        fwrite(fid, 0, 'int');        % bonf
        fwrite(fid, [255,0,0], 'uchar');   % col1
        fwrite(fid, [255,255,0], 'uchar'); % col2
        fwrite(fid, 1, 'char');       % enable SMP
        fwrite(fid, 1, 'float');      % transparency
        fwrite(fid, 0, 'char');       % name

        fwrite(fid, siz, 'int');      % original size
        fwrite(fid, minx-1, 'int');
        fwrite(fid, maxx-1, 'int');
        fwrite(fid, miny-1, 'int');
        fwrite(fid, maxy-1, 'int');
        fwrite(fid, minz-1, 'int');
        fwrite(fid, maxz-1, 'int');
        fwrite(fid, 1, 'int');        % resolution

        % write the data
        fwrite(fid, data(minx:maxx,miny:maxy,minz:maxz), 'float');
    end
    fclose(fid);

  case 'vmr'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write in BrainVoyager VMR format
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen(sprintf('%s.vmr',cfg.filename),'w');
    if fid < 0,
      error('Cannot write to file %s.vmr\n',cfg.filename);
    end

    % data should be scaled between 0 and 225
    data = data - min(data(:));
    data = round(225*data./max(data(:)));

    % write the header
    fwrite(fid, siz, 'ushort');
    % write the data
    fwrite(fid, data, 'uint8');
    fclose(fid);
  case 'analyze'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write in Analyze format, using some functions from Darren Webbers toolbox
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    avw = avw_hdr_make;

    % specify the image data and dimensions
    avw.hdr.dime.dim(2:4) = siz;
    avw.img = data;

    % orientation 0 means transverse unflipped (axial, radiological)
    % X direction first,  ft_progressing from patient right to left,
    % Y direction second, ft_progressing from patient posterior to anterior,
    % Z direction third,  ft_progressing from patient inferior to superior.
    avw.hdr.hist.orient = 0;

    % specify voxel size
    avw.hdr.dime.pixdim(2:4) = [1 1 1];
    % FIXME, this currently does not work due to all flipping and permuting
    % resx = x(2)-x(1);
    % resy = y(2)-y(1);
    % resz = z(2)-z(1);
    % avw.hdr.dime.pixdim(2:4) = [resy resx resz];

    % specify the data type
    switch lower(cfg.datatype)
      case 'bit1'
        avw.hdr.dime.datatype = 1;
        avw.hdr.dime.bitpix   = 1;
      case 'uint8'
        avw.hdr.dime.datatype = 2;
        avw.hdr.dime.bitpix   = 8;
      case 'int16'
        avw.hdr.dime.datatype = 4;
        avw.hdr.dime.bitpix   = 16;
      case 'int32'
        avw.hdr.dime.datatype = 8;
        avw.hdr.dime.bitpix   = 32;
      case 'float'
        avw.hdr.dime.datatype = 16;
        avw.hdr.dime.bitpix   = 32;
      case 'double'
        avw.hdr.dime.datatype = 64;
        avw.hdr.dime.bitpix   = 64;
      otherwise
        error('unknown datatype');
    end

    % write the header and image data
    avw_img_write(avw, cfg.filename, [], 'ieee-le');

  case 'nifti'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write in nifti format, using functions from  the SPM8 toolbox
    % this format supports a homogenous transformation matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pathstr, name, ext] = fileparts(cfg.filename);
    if isempty(ext)
      cfg.filename = [cfg.filename,'.nii'];
    end
    ft_write_mri(cfg.filename, data, 'dataformat', 'nifti', 'transform', transform, 'spmversion', 'SPM8');

  case 'nifti_img'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write in nifti dual file format, using functions from  the SPM8 toolbox
    % this format supports a homogenous transformation matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pathstr, name, ext] = fileparts(cfg.filename);
    if isempty(ext)
      cfg.filename = [cfg.filename,'.img'];
    end
    ft_write_mri(cfg.filename, data, 'dataformat', 'nifti', 'transform', transform, 'spmversion', 'SPM8');

  case 'analyze_spm'
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write in analyze format, using functions from  the SPM8 toolbox
    % this format supports a homogenous transformation matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pathstr, name, ext] = fileparts(cfg.filename);
    if isempty(ext)
      cfg.filename = [cfg.filename,'.img'];
    end
    ft_write_mri(cfg.filename, data, 'dataformat', 'analyze', 'transform', transform, 'spmversion', 'SPM2');

  case {'mgz' 'mgh'}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % write in freesurfer_mgz format, using functions from  the freesurfer toolbox
    % this format supports a homogenous transformation matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ispc && strcmp(cfg.filetype, 'mgz')
      warning('Saving in .mgz format is not possible on a PC, saving in .mgh format instead');
      cfg.filetype = 'mgh';
    end
    [pathstr, name, ext] = fileparts(cfg.filename);
    if isempty(ext)
      cfg.filename = [cfg.filename,'.',cfg.filetype];
    end
    ft_write_mri(cfg.filename, data, 'dataformat', cfg.filetype, 'transform', transform);


  otherwise
    fprintf('unknown fileformat\n');
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble previous volume
ft_postamble provenance

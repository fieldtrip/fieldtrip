function [mri] = ft_read_mri(filename, varargin)

% FT_READ_MRI reads anatomical and functional MRI data from different file formats.
% The output data is structured in such a way that it is compatible with
% FT_DATATYPE_VOLUME.
%
% Use as
%   [mri] = ft_read_mri(filename)
%
% Additional options should be specified in key-value pairs and can be
%   'dataformat' = string specifying the file format, determining the low-
%                  level reading routine to be used. If no explicit format
%                  is given, it is determined automatically from the filename.
%   'outputfield' = string specifying the name of the field in the
%                  structure in which the numeric data is stored. The
%                  default is 'anatomy'
%
% The following values apply for the dataformat
%   'afni_head'/'afni_brik'      uses AFNI code
%   'analyze_img'/'analyze_hdr'  uses SPM code
%   'analyze_old'                uses Darren Webber's code
%   'asa_mri'
%   'ctf_mri'
%   'ctf_mri4'
%   'ctf_svl'
%   'dicom'                      uses FreeSurfer code
%   'dicom_old'                  uses own code
%   'freesurfer_mgh'             uses FreeSurfer code
%   'freesurfer_mgz'             uses FreeSurfer code
%   'matlab'                     assumes a MATLAB *.mat file containing a mri structure according to FT_DATATYPE_VOLUME
%   'minc'                       uses SPM (<= version SPM5)
%   'neuromag_fif'               uses MNE toolbox
%   'neuromag_fif_old'           uses meg-pd toolbox
%   'nifti'                      uses FreeSurfer code
%   'nifti_fsl'                  uses FreeSurfer code
%   'nifti_spm'                  uses SPM
%   'yokogawa_mri'
%   'mrtrix_mif'                 uses mrtrix code
%
% The following MRI file formats are supported
%   CTF (*.svl, *.mri version 4 and 5)
%   NIFTi (*.nii) and zipped NIFTi (*.nii.gz)
%   Analyze (*.img, *.hdr)
%   DICOM (*.dcm, *.ima)
%   AFNI (*.head, *.brik)
%   FreeSurfer (*.mgz, *.mgh)
%   MINC (*.mnc)
%   Neuromag/Elekta/Megin (*.fif)
%   ANT - Advanced Neuro Technology (*.mri)
%   Yokogawa (*.mrk, incomplete)
%   Mrtrix image format (*.mif)
%
% If you have a series of DICOM files, please provide the name of any of the files
% in the series (e.g. the first one). The other files will be found automatically.
%
% The output MRI may have a homogenous transformation matrix that converts
% the coordinates of each voxel (in xgrid/ygrid/zgrid) into head
% coordinates.
%
% See also FT_DATATYPE_VOLUME, FT_WRITE_MRI, FT_READ_DATA, FT_READ_HEADER, FT_READ_EVENT

% Undocumented options
%   'fixel2voxel' = string, operation to apply to the fixels belonging to  the same voxel (only for *.mif). 'max' (default), 'min', 'mean'
%   'indexfile'   = string, pointing to a fixel index file, if not present in the same directory as the functional data

% Copyright (C) 2008-2020, Robert Oostenveld & Jan-Mathijs Schoffelen
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

% optionally get the data from the URL and make a temporary local copy
filename = fetch_url(filename);

% get the options
dataformat  = ft_getopt(varargin, 'dataformat');
outputfield = ft_getopt(varargin, 'outputfield', 'anatomy');

% the following is added for backward compatibility of using 'format' rather than 'dataformat'
format    = ft_getopt(varargin, 'format');
if ~isempty(format)
  ft_warning('the option ''format'' will be deprecated soon, please use ''dataformat'' instead');
  if isempty(dataformat)
    dataformat  = format;
  end
end

if isempty(dataformat)
  % only do the autodetection if the format was not specified
  dataformat = ft_filetype(filename);
end

if strcmp(dataformat, 'compressed') || (strcmp(dataformat, 'freesurfer_mgz') && ispc)
  % the file is compressed, unzip on the fly, freesurfer mgz files get
  % special treatment on a pc
  inflated = true;
  filename = inflate_file(filename);
  if strcmp(dataformat, 'freesurfer_mgz')
    filename_old = filename;
    filename     = [filename '.mgh'];
    movefile(filename_old, filename);
  end
  dataformat = ft_filetype(filename);
else
  inflated = false;
end

% test whether the file exists
if ~exist(filename, 'file')
  ft_error('file ''%s'' does not exist', filename);
end

% test for the presence of some external functions from other toolboxes
hasspm2  = ft_hastoolbox('spm2');    % see http://www.fil.ion.ucl.ac.uk/spm/
hasspm5  = ft_hastoolbox('spm5');    % see http://www.fil.ion.ucl.ac.uk/spm/
hasspm8  = ft_hastoolbox('spm8');    % see http://www.fil.ion.ucl.ac.uk/spm/
hasspm12 = ft_hastoolbox('spm12');   % see http://www.fil.ion.ucl.ac.uk/spm/

switch dataformat
  case 'ctf_mri'
    [img, hdr] = read_ctf_mri(filename);
    transform = hdr.transformMRI2Head;
    coordsys  = 'ctf';

  case 'ctf_mri4'
    [img, hdr] = read_ctf_mri4(filename);
    transform = hdr.transformMRI2Head;
    coordsys  = 'ctf';

  case 'ctf_svl'
    [img, hdr] = read_ctf_svl(filename);
    transform = hdr.transform;

  case 'asa_mri'
    [img, seg, hdr] = read_asa_mri(filename);
    transform = hdr.transformMRI2Head;
    if isempty(seg)
      % in case seg exists it will be added to the output
      clear seg
    end

  case 'minc'
    if ~(hasspm2 || hasspm5)
      fprintf('the SPM2 or SPM5 toolbox is required to read *.mnc files\n');
      ft_hastoolbox('spm2', 1);
    end
    % use the functions from SPM
    hdr = spm_vol_minc(filename);
    img = spm_read_vols(hdr);
    transform = hdr.mat;

  case 'nifti_spm'
    if ~(hasspm5 || hasspm8 || hasspm12)
      fprintf('the SPM5 or newer toolbox is required to read *.nii files\n');
      ft_hastoolbox('spm12', 1);
    end
    % use the functions from SPM
    hdr = spm_vol_nifti(filename);
    img = double(hdr.private.dat);
    %img = spm_read_vols(hdr);
    transform = hdr.mat;

  case {'analyze_img' 'analyze_hdr'}
    if ~(hasspm8 || hasspm12)
      fprintf('the SPM8 or newer toolbox is required to read analyze files\n');
      ft_hastoolbox('spm8up', 1);
    end

    % use the image file instead of the header
    filename((end-2):end) = 'img';
    % use the functions from SPM to read the Analyze MRI
    hdr = spm_vol(filename);
    img = spm_read_vols(hdr);
    transform = hdr.mat;

  case 'analyze_old'
    % use the functions from Darren Weber's mri_toolbox to read the Analyze MRI
    ft_hastoolbox('mri', 1);     % from Darren Weber, see http://eeg.sourceforge.net/

    avw = avw_img_read(filename, 0); % returned volume is LAS*
    img = avw.img;
    hdr = avw.hdr;
    % The default Analyze orientation is axial unflipped (LAS*), which means
    % that the resulting volume is according to the radiological convention.
    % Most other fMRI and EEG/MEG software (except Mayo/Analyze) uses
    % neurological conventions and a right-handed coordinate system, hence
    % the first axis of the 3D volume (right-left) should be flipped to make
    % the coordinate system comparable to SPM
    ft_warning('flipping 1st dimension (L-R) to obtain volume in neurological convention');
    img = flipdim(img, 1);

    transform      = diag(hdr.dime.pixdim(2:4));
    transform(4,4) = 1;

  case {'afni_brik' 'afni_head'}
    % needs afni
    ft_hastoolbox('afni', 1);    % see http://afni.nimh.nih.gov/

    [err, img, hdr, ErrMessage] = BrikLoad(filename);
    if err
      ft_error('could not read AFNI file');
    end

    % FIXME: this should be checked, but I only have a single BRIK file
    % construct the homogenous transformation matrix that defines the axes
    ft_warning('homogenous transformation might be incorrect for AFNI file');
    transform        = eye(4);
    transform(1:3,4) = hdr.ORIGIN(:);
    transform(1,1)   = hdr.DELTA(1);
    transform(2,2)   = hdr.DELTA(2);
    transform(3,3)   = hdr.DELTA(3);

    % FIXME: I am not sure about the "RAI" image orientation
    img = flipdim(img,1);
    img = flipdim(img,2);
    dim = size(img);
    transform(1,4) = -dim(1) - transform(1,4);
    transform(2,4) = -dim(2) - transform(2,4);

  case 'neuromag_fif'
    % needs mne toolbox
    ft_hastoolbox('mne', 1);

    % use the mne functions to read the Neuromag MRI
    hdr = fiff_read_mri(filename);
    img_t = cat(3, hdr.slices.data);
    img = permute(img_t,[2 1 3]);
    hdr.slices = rmfield(hdr.slices, 'data'); % remove the image data to save memory

    % information below is from MNE - fiff_define_constants.m
    % coordinate system 4 - is the MEG head coordinate system (fiducials)
    % coordinate system 5 - is the MRI coordinate system
    % coordinate system 2001 - MRI voxel coordinates
    % coordinate system 2002 - Surface RAS coordinates (is mainly vertical
    %                                     shift, no rotation to 2001)
    % MEG sensor positions come in system 4
    % MRI comes in system 2001

    transform = eye(4);
    if isfield(hdr, 'trans') && issubfield(hdr.trans, 'trans')
      if (hdr.trans.from == 4) && (hdr.trans.to == 5)
        transform = hdr.trans.trans;
      else
        ft_warning('W: trans does not transform from 4 to 5.');
        ft_warning('W: Please check the MRI fif-file');
      end
    else
      ft_warning('W: trans structure is not defined.');
      ft_warning('W: Maybe coregistration is missing?');
    end
    if isfield(hdr, 'voxel_trans') && issubfield(hdr.voxel_trans, 'trans')
      % centers the coordinate system
      % and switches from mm to m
      if (hdr.voxel_trans.from == 2001) && (hdr.voxel_trans.to == 5)
        % matlab_shift compensates for the different index conventions
        % between C and matlab

        % the lines below is old code (prior to Jan 3, 2013) and only works with
        % 1 mm resolution MRIs
        %matlab_shift = [ 0 0 0 0.001; 0 0 0 -0.001; 0 0 0 0.001; 0 0 0 0];
        % transform transforms from 2001 to 5 and further to 4
        %transform = transform\(hdr.voxel_trans.trans+matlab_shift);

        % the lines below should work with arbitrary resolution
        matlab_shift = eye(4);
        matlab_shift(1:3,4) = [-1,-1,-1];
        transform = transform\(hdr.voxel_trans.trans * matlab_shift);

        coordsys  = 'neuromag';
        mri.unit  = 'm';
      else
        ft_warning('W: voxel_trans does not transform from 2001 to 5.');
        ft_warning('W: Please check the MRI fif-file');
      end
    else
      ft_warning('W: voxel_trans structure is not defined.');
      ft_warning('W: Please check the MRI fif-file');
    end

  case 'neuromag_fif_old'
    % needs meg_pd toolbox
    ft_hastoolbox('meg-pd', 1);

    % use the meg_pd functions to read the Neuromag MRI
    [img,coords] = loadmri(filename);
    dev = loadtrans(filename,'MRI','HEAD');
    transform  = dev*coords;
    hdr.coords = coords;
    hdr.dev    = dev;

  case 'dicom'
    % this seems to return a right-handed volume with the transformation
    % matrix stored in the file headers.

    % needs the freesurfer toolbox
    ft_hastoolbox('freesurfer', 1);
    [dcmdir,junk1,junk2] = fileparts(filename);
    if isempty(dcmdir),
      dcmdir = '.';
    end
    [img,transform,hdr,mr_params] = load_dicom_series(dcmdir,dcmdir,filename);
    transform = vox2ras_0to1(transform);

  case 'dicom_old'
    % this does not necessarily return a right-handed volume and only a
    % transformation-matrix with the voxel size

    % this uses the Image processing toolbox
    % the DICOM file probably represents a stack of slices, possibly even multiple volumes
    orig = dicominfo(filename);
    dim(1) = orig.Rows;
    dim(2) = orig.Columns;

    [p, f] = fileparts(filename);

    % this works for the Siemens scanners at the FCDC
    tok = tokenize(f, '.');
    for i=5:length(tok)
      tok{i} = '*';
    end
    filename = sprintf('%s.', tok{:});  % reconstruct the filename with wildcards and '.' between the segments
    filename = filename(1:end-1);       % remove the last '.'
    dirlist  = dir(fullfile(p, filename));
    dirlist  = {dirlist.name};

    if isempty(dirlist)
      % this is for the Philips data acquired at KI
      ft_warning('could not determine list of dicom files, trying with *.dcm');
      dirlist  = dir(fullfile(p, '*.dcm'));
      dirlist  = {dirlist.name};
    end

    if isempty(dirlist)
      ft_warning('could not determine list of dicom files, trying with *.ima');
      dirlist  = dir(fullfile(p, '*.ima'));
      dirlist  = {dirlist.name};
    end

    if length(dirlist)==1
      % try something else to get a list of all the slices
      dirlist = dir(fullfile(p, '*'));
      dirlist = {dirlist(~[dirlist.isdir]).name};
    end

    keep = false(1, length(dirlist));
    for i=1:length(dirlist)
      filename = char(fullfile(p, dirlist{i}));
      if ~strcmp(dataformat, 'dicom_old')
        keep(i) = false;
        fprintf('skipping ''%s'' because of incorrect filetype\n', filename);
      end
      % read the header information
      info     = dicominfo(filename);
      if info.SeriesNumber~=orig.SeriesNumber
        keep(i) = false;
        fprintf('skipping ''%s'' because of different SeriesNumber\n', filename);
      else
        keep(i) = true;
        hdr(i)  = info;
      end
    end
    % remove the files that were skipped
    hdr     = hdr(keep);
    dirlist = dirlist(keep);

    % pre-allocate enough space for the subsequent slices
    dim(3) = length(dirlist);
    img    = zeros(dim(1), dim(2), dim(3));
    for i=1:length(dirlist)
      filename = char(fullfile(p, dirlist{i}));
      fprintf('reading image data from ''%s''\n', filename);
      img(:,:,i) = dicomread(hdr(i));
    end

    % reorder the slices
    [z, indx]   = sort(cell2mat({hdr.SliceLocation}));
    hdr = hdr(indx);
    img = img(:,:,indx);

    try
      % construct a homgeneous transformation matrix that performs the scaling from voxels to mm
      dx = hdr(1).PixelSpacing(1);
      dy = hdr(1).PixelSpacing(2);
      dz = hdr(2).SliceLocation - hdr(1).SliceLocation;
      transform = eye(4);
      transform(1,1) = dx;
      transform(2,2) = dy;
      transform(3,3) = dz;
    end

  case {'nifti', 'freesurfer_mgz', 'freesurfer_mgh', 'nifti_fsl'}

    ft_hastoolbox('freesurfer', 1);
    tmp = MRIread(filename);
    ndims = numel(size(tmp.vol));
    if ndims==3
      img = permute(tmp.vol, [2 1 3]); %FIXME although this is probably correct
      %see the help of MRIread, anecdotally columns and rows seem to need a swap
      %in order to match the transform matrix (alternatively a row switch of the
      %latter can be done)
    elseif ndims==4
      img = permute(tmp.vol, [2 1 3 4]);
    end
    hdr = rmfield(tmp, 'vol');
    transform = tmp.vox2ras1;

  case 'yokogawa_mri'
    ft_hastoolbox('yokogawa', 1);
    fid = fopen_or_error(filename, 'rb');
    mri_info = GetMeg160MriInfoM(fid);
    patient_info = GetMeg160PatientInfoFromMriFileM(fid);
    [data_style, model, marker, image_parameter, normalize, besa_fiducial_point] = GetMeg160MriFileHeaderInfoM(fid);
    fclose(fid);

    % gather all meta-information
    hdr.mri_info = mri_info;
    hdr.patient_info = patient_info;
    hdr.data_style = data_style;
    hdr.model = model;
    hdr.marker = marker;
    hdr.image_parameter = image_parameter;
    hdr.normalize = normalize;
    hdr.besa_fiducial_point = besa_fiducial_point;

    ft_error('FIXME yokogawa_mri implementation is incomplete');

  case 'matlab'
    mri = loadvar(filename, 'mri');

  case {'mif' 'mrtrix_mif'}
    ft_hastoolbox('mrtrix', 1);
    
    tmp = read_mrtrix(filename);
    
    % check if it's sparse fixeldata
    isfixel = numel(tmp.dim==3)&&tmp.dim(3)==1;
      
    if ~isfixel
      mri.hdr     = removefields(tmp, {'data'});
      mri.(outputfield) = tmp.data;
      mri.dim     = tmp.dim(1:length(size(tmp.data)));
      mri.transform = tmp.transform;
      mri.transform(1:3,1:3) = diag(tmp.vox(1:3))*mri.transform(1:3,1:3);
    else
      fix2vox_fun = ft_getopt(varargin, 'fixel2voxel', 'max');
      indexfile   = ft_getopt(varargin, 'indexfile');
      if isempty(indexfile)
        % assume the index file to be in the same directory as the data file
        [p,f,e]   = fileparts(filename);
        indexfile = fullfile(p,'index.mif');
      end
      index     = read_mrtrix(indexfile);
      tmpdata   = reshape(index.data, [], 2);
      
      vox_index = find(tmpdata(:,1)>0);
      num_index = tmpdata(vox_index,1);
      fix_index = tmpdata(vox_index,2)+1;
      
      % create a mapping matrix of fixel2voxel -> currently this only works
      % for scalar fixel data.
      tmpdata = zeros(numel(num_index), max(num_index));
      for k = 1:numel(num_index)
        tmpdata(k, 1:num_index(k)) = fix_index(k)-1 + (1:num_index(k));
      end
      tmpdata(tmpdata==0) = nan;
      tmpdata             = tmpdata.'; % transpose is intended
      tmpdata(isfinite(tmpdata)) = tmp.data(tmpdata(isfinite(tmpdata)));
      
      switch fix2vox_fun
        case 'magmax'
          tmpx    = nanmin(tmpdata,[],1).';
          tmpdata = nanmax(tmpdata,[],1).';
          tmpdata(abs(tmpx)>tmpdata) = tmpx(abs(tmpx)>tmpdata);
        case 'max'
          tmpdata = nanmax(tmpdata,[],1).';
        case 'min'
          tmpdata = nanmin(tmpdata,[],1).';
        case 'mean'
          tmpdata = nanmean(tmpdata,1).';
        case 'none'
          tmpdata = tmpdata.';
        otherwise
          ft_error('unsupported fixel2voxel operation requested');
      end
      
      mri.hdr  = removefields(tmp, {'data'});
      mri.(outputfield) = zeros([prod(index.dim(1:3)) size(tmpdata,2)]);
      mri.(outputfield)(vox_index,:) = tmpdata;
      mri.(outputfield)   = reshape(mri.(outputfield), [index.dim(1:3) size(tmpdata,2)]);
      mri.dim       = size(mri.(outputfield)) ;%index.dim(1:length(size(tmp.data)));
      mri.transform = tmp.transform;
      mri.transform(1:3,1:3) = diag(tmp.vox(1:3))*mri.transform(1:3,1:3);
    end
    
  otherwise
    ft_error('unrecognized filetype ''%s'' for ''%s''', dataformat, filename);
end

if exist('img', 'var')
  % set up the axes of the volume in voxel coordinates
  %nx = size(img,1);
  %ny = size(img,2);
  %nz = size(img,3);
  mri.dim = size(img); %[nx ny nz];
  % store the anatomical data
  mri.(outputfield) = img;
end

if exist('seg', 'var')
  % store the segmented data
  mri.seg = seg;
end

if exist('hdr', 'var')
  % store the header with all file format specific details
  mri.hdr = hdr;
end

try
  % store the homogenous transformation matrix if present
  mri.transform = transform;
end

try
  % try to determine the units of the coordinate system
  mri = ft_determine_units(mri);
end

try
  % try to add a descriptive label for the coordinate system
  mri.coordsys = coordsys;
end

if inflated
  % compressed file has been unzipped on the fly, clean up
  delete(filename);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function value = loadvar(filename, varname)
var = whos('-file', filename);
if length(var)==1
  filecontent = load(filename); % read the one variable in the file, regardless of how it is called
  value       = filecontent.(var.name);
  clear filecontent
else
  filecontent = load(filename, varname);
  value       = filecontent.(varname);  % read the variable named according to the input specification
  clear filecontent
end

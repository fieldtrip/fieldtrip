function [mri] = ft_read_mri(filename, varargin)

% FT_READ_MRI reads anatomical and functional MRI data from different
% file formats. The output data is structured in such a way that it is
% comparable to a FieldTrip source reconstruction.
%
% Use as
%   [mri] = ft_read_mri(filename)
%
% The output MRI may have a homogenous transformation matrix that converts
% the coordinates of each voxel (in xgrid/ygrid/zgrid) into head
% coordinates.
%
% The following MRI file formats are supported
%   CTF - VSM MedTech (*.svl, *.mri version 4 and 5)
%   NIFTi (*.nii)
%   Analyze (*.img, *.hdr)
%   DICOM (*.dcm, *.ima)
%   AFNI (*.head, *.brik)
%   FreeSurfer (*.mgz, *.mgh)
%   MINC (*.mnc)
%   Neuromag - Elekta (*.fif)
%   ANT - Advanced Neuro Technology (*.mri)
%   Yokogawa (*.mrk, incomplete)
%
% See also FT_WRITE_MRI, FT_READ_DATA, FT_READ_HEADER, FT_READ_EVENT

% Copyright (C) 2008-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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
mriformat = ft_getopt(varargin, 'format'); % FIXME this is inconsistent with ft_read_mri, which uses 'dataformat'

if isempty(mriformat)
  % only do the autodetection if the format was not specified
  mriformat = ft_filetype(filename);
end

% test whether the file exists
if ~exist(filename, 'file')
  error('file ''%s'' does not exist', filename);
end

% test for the presence of some external functions from other toolboxes
hasmri   = ft_hastoolbox('mri');     % from Darren Weber, see http://eeg.sourceforge.net/
hasspm2  = ft_hastoolbox('spm2');    % see http://www.fil.ion.ucl.ac.uk/spm/
hasspm5  = ft_hastoolbox('spm5');    % see http://www.fil.ion.ucl.ac.uk/spm/
hasspm8  = ft_hastoolbox('spm8');    % see http://www.fil.ion.ucl.ac.uk/spm/
hasspm12 = ft_hastoolbox('spm12');   % see http://www.fil.ion.ucl.ac.uk/spm/
hasspm   = (hasspm2 || hasspm5 || hasspm8 || hasspm12);
hasafni  = ft_hastoolbox('afni');    % see http://afni.nimh.nih.gov/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mriformat, 'ctf_mri')
  [img, hdr] = read_ctf_mri(filename);
  transform = hdr.transformMRI2Head;
  coordsys  = 'ctf';
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mriformat, 'ctf_mri4')
  [img, hdr] = read_ctf_mri4(filename);
  transform = hdr.transformMRI2Head;
  coordsys  = 'ctf';
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mriformat, 'ctf_svl')
  [img, hdr] = read_ctf_svl(filename);
  transform = hdr.transform;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mriformat, 'asa_mri')
  [img, seg, hdr] = read_asa_mri(filename);
  transform = hdr.transformMRI2Head;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mriformat, 'minc')
  if ~(hasspm2 || hasspm5)
    fprintf('the SPM2 or SPM5 toolbox is required to read *.mnc files\n');
    ft_hastoolbox('spm2',1);
  end
  % use the functions from SPM
  hdr = spm_vol_minc(filename);
  img = spm_read_vols(hdr);
  transform = hdr.mat;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % USE FREESURFER CODE FOR THE READING OF NIFTI-FILES: THAT CODE ALSO
  % DEALS WITH 4D NIFTIs
elseif strcmp(mriformat, 'nifti_spm')
  if ~(hasspm5 || hasspm8 || hasspm12)
    fprintf('the SPM5 or newer toolbox is required to read *.nii files\n');
    ft_hastoolbox('spm8',1);
  end
  % use the functions from SPM
  hdr = spm_vol_nifti(filename);
  img = spm_read_vols(hdr);
  transform = hdr.mat;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (strcmp(mriformat, 'analyze_img') || strcmp(mriformat, 'analyze_hdr')) && hasspm
  % use the image file instead of the header
  filename((end-2):end) = 'img';
  % use the functions from SPM to read the Analyze MRI
  hdr = spm_vol(filename);
  img = spm_read_vols(hdr);
  transform = hdr.mat;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (strcmp(mriformat, 'analyze_hdr') || strcmp(mriformat, 'analyze_img')) && hasmri
  % use the functions from Darren Weber's mri_toolbox to read the Analyze MRI
  avw = avw_img_read(filename, 0); % returned volume is LAS*
  img = avw.img;
  hdr = avw.hdr;
  % The default Analyze orientation is axial unflipped (LAS*), which means
  % that the resulting volume is according to the radiological convention.
  % Most other fMRI and EEG/MEG software (except Mayo/Analyze) uses
  % neurological conventions and a right-handed coordinate system, hence
  % the first axis of the 3D volume (right-left) should be flipped to make
  % the coordinate system comparable to SPM
  warning('flipping 1st dimension (L-R) to obtain volume in neurological convention');
  img = flipdim(img, 1);
  
  transform      = diag(hdr.dime.pixdim(2:4));
  transform(4,4) = 1;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (strcmp(mriformat, 'afni_brik') || strcmp(mriformat, 'afni_head')) && hasafni
  [err, img, hdr, ErrMessage] = BrikLoad(filename);
  if err
    error('could not read AFNI file');
  end
  
  % FIXME: this should be checked, but I only have a single BRIK file
  % construct the homogenous transformation matrix that defines the axes
  warning('homogenous transformation might be incorrect for AFNI file');
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
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mriformat, 'neuromag_fif') && ft_hastoolbox('mne')
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
      warning('W: trans does not transform from 4 to 5.');
      warning('W: Please check the MRI fif-file');
    end
  else
    warning('W: trans structure is not defined.');
    warning('W: Maybe coregistration is missing?');
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
      warning('W: voxel_trans does not transform from 2001 to 5.');
      warning('W: Please check the MRI fif-file');
    end
  else
    warning('W: voxel_trans structure is not defined.');
    warning('W: Please check the MRI fif-file');
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mriformat, 'neuromag_fif') && ft_hastoolbox('meg_pd')
  % use the meg_pd functions to read the Neuromag MRI
  [img,coords] = loadmri(filename);
  dev = loadtrans(filename,'MRI','HEAD');
  transform  = dev*coords;
  hdr.coords = coords;
  hdr.dev    = dev;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mriformat, 'neuromag_fif')
  error('reading MRI data from a fif file requires either the MNE toolbox or the meg_pd toolbox to be installed');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mriformat, 'dicom')
  % this seems to return a right-handed volume with the transformation
  % matrix stored in the file headers.
  
  % use the freesurfer toolbox
  ft_hastoolbox('freesurfer', 1);
  [dcmdir,junk1,junk2] = fileparts(filename);
  if isempty(dcmdir),
    dcmdir = '.';
  end
  [img,transform,hdr,mr_params] = load_dicom_series(dcmdir,dcmdir,filename);
  transform = vox2ras_0to1(transform);
  
elseif strcmp(mriformat, 'dicom_old')
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
  
  if length(dirlist)==1
    % try something else to get a list of all the slices
    dirlist = dir(fullfile(p, '*'));
    dirlist = {dirlist(~[dirlist.isdir]).name};
  end
  
  keep = false(1, length(dirlist));
  for i=1:length(dirlist)
    filename = char(fullfile(p, dirlist{i}));
    if ~strcmp(mriformat, 'dicom')
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
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif any(strcmp(mriformat, {'nifti', 'freesurfer_mgz', 'freesurfer_mgh', 'nifti_fsl'}))
  if strcmp(mriformat, 'freesurfer_mgz') && ispc
    error('Compressed .mgz files cannot be read on a PC');
  end
  
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
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mriformat, 'yokogawa_mri')
  ft_hastoolbox('yokogawa', 1);
  fid = fopen(filename, 'rb');
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
  
  error('FIXME yokogawa_mri implementation is incomplete');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(mriformat, 'matlab')
  mri = loadvar(filename, 'mri');
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  error(sprintf('unrecognized filetype ''%s'' for ''%s''', mriformat, filename));
end

if exist('img', 'var')
  % set up the axes of the volume in voxel coordinates
  nx = size(img,1);
  ny = size(img,2);
  nz = size(img,3);
  mri.dim = [nx ny nz];
  % store the anatomical data
  mri.anatomy = img;
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
  mri = ft_convert_units(mri);
end

try
  % try to add a descriptive label for the coordinate system
  mri.coordsys = coordsys;
end

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

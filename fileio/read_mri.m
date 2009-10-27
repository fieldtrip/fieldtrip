function [mri] = read_mri(filename)

% READ_MRI reads anatomical and functional MRI data from different
% file formats. The output data is structured in such a way that it is
% comparable to a FieldTrip source reconstruction.
%
% Use as
%   [mri] = read_mri(filename)
%
% The output MRI may have a homogenous transformation matrix that converts
% the coordinates of each voxel (in xgrid/ygrid/zgrid) into head
% coordinates.
%
% See also READ_DATA, READ_HEADER, READ_EVENT

% Copyright (C) 2004-2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fieldtripdefs

% test for the presence of some external functions from other toolboxes
hasmri  = hastoolbox('mri');     % from Darren Weber, see http://eeg.sourceforge.net/
hasspm2 = hastoolbox('spm2');    % see http://www.fil.ion.ucl.ac.uk/spm/
hasspm5 = hastoolbox('spm5');    % see http://www.fil.ion.ucl.ac.uk/spm/
hasspm8 = hastoolbox('spm8');    % see http://www.fil.ion.ucl.ac.uk/spm/
hasspm = (hasspm2 || hasspm5 || hasspm8);
hasafni = hastoolbox('afni');    % see http://afni.nimh.nih.gov/

% test whether the file exists
if ~exist(filename)
  error(sprintf('file ''%s'' does not exist', filename));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if filetype(filename, 'ctf_mri')
  [img, hdr] = read_ctf_mri(filename);
  transform = hdr.transformMRI2Head;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filetype(filename, 'ctf_mri4')
  [img, hdr] = read_ctf_mri4(filename);
  transform = hdr.transformMRI2Head;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filetype(filename, 'asa_mri')
  [img, seg, hdr] = read_asa_mri(filename);
  transform = hdr.transformMRI2Head;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filetype(filename, 'minc')
  if ~hasspm
    error('the SPM2 or SPM5 toolbox is required to read *.mnc files');
  end
  % use the functions from SPM
  hdr = spm_vol_minc(filename);
  img = spm_read_vols(hdr);
  transform = hdr.mat;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filetype(filename, 'nifti')
  if ~(hasspm5 || hasspm8)
    error('the SPM5 or SPM8 toolbox is required to read *.nii files');
  end
  % use the functions from SPM
  hdr = spm_vol_nifti(filename);
  img = spm_read_vols(hdr);
  transform = hdr.mat;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (filetype(filename, 'analyze_img') || filetype(filename, 'analyze_hdr')) && hasspm
  % use the image file instead of the header
  filename((end-2):end) = 'img';
  % use the functions from SPM to read the Analyze MRI
  hdr = spm_vol(filename);
  img = spm_read_vols(hdr);
  transform = hdr.mat;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (filetype(filename, 'analyze_hdr') || filetype(filename, 'analyze_img')) && hasmri
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
  % FIXME: here I should also implement a homogenous transformation matrix,
  % using the voxel dimensions that are specified in hdr.dime.pixdim

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (filetype(filename, 'afni_brik') || filetype(filename, 'afni_head')) && hasafni
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
elseif filetype(filename, 'neuromag_fif') && hastoolbox('mne')
  % use the mne functions to read the Neuromag MRI
  hdr = fiff_read_mri(filename);
  img = cat(3, hdr.slices.data);
  hdr.slices = rmfield(hdr.slices, 'data'); % remove the image data to save memory
  % hmm, which transformation matrix should I use?
  if issubfield(hdr.voxel_trans, 'trans')
    transform = hdr.voxel_trans.trans;
  elseif issubfield(hdr.trans, 'trans')
    transform = hdr.trans.trans;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filetype(filename, 'neuromag_fif') && hastoolbox('meg_pd')
  % use the meg_pd functions to read the Neuromag MRI
  [img,coords] = loadmri(filename);
  dev = loadtrans(filename,'MRI','HEAD');
  transform  = dev*coords;
  hdr.coords = coords;
  hdr.dev    = dev;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filetype(filename, 'neuromag_fif')
  error('reading MRI data from a fif file requires either the MNE toolbox or the meg_pd toolbox to be installed');

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif filetype(filename, 'dicom')
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
    if ~filetype(filename, 'dicom')
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
    % construct a homgenous transformation matrix that performs the scaling from voxels to mm
    dx = hdr(1).PixelSpacing(1);
    dy = hdr(1).PixelSpacing(2);
    dz = hdr(2).SliceLocation - hdr(1).SliceLocation;
    transform = eye(4);
    transform(1,1) = dx;
    transform(2,2) = dy;
    transform(3,3) = dz;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  error(sprintf('unrecognized filetype of ''%s''', filename));
end

% set up the axes of the volume in voxel coordinates
nx = size(img,1);
ny = size(img,2);
nz = size(img,3);
mri.dim = [nx ny nz];
% store the anatomical data
mri.anatomy = img;
% store the header with all fileformat specific details
mri.hdr = hdr;
try
  % if present, store the homogenous transformation matrix
  mri.transform = transform;
end


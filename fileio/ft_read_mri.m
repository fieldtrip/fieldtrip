function [mri] = ft_read_mri(filename, varargin)

% FT_READ_MRI reads anatomical and functional MRI data from different file formats.
% The output data is structured in such a way that it is compatible with
% FT_DATATYPE_VOLUME.
%
% Use as
%   [mri] = ft_read_mri(filename, ...)
%
% Additional options should be specified in key-value pairs and can include
%   'dataformat'  = string specifying the file format, determining the low-level
%                   reading routine to be used. If no explicit format is given,
%                   it is determined automatically from the filename.
%   'volumes'     = vector with the volume indices to read from a 4D nifti (only for 'nifti_spm')
%   'outputfield' = string specifying the name of the field in the structure in which the
%                   numeric data is stored (only for 'mrtrix_mif', default = 'anatomy')
%   'fixel2voxel' = string, the operation to apply to the fixels belonging to the
%                  same voxel, can be 'max', 'min', 'mean' (only for 'mrtrix_mif', default = 'max')
%   'indexfile'   = string, pointing to a fixel index file, if not present in the same directory
%                   as the functional data (only for 'mrtrix_mif')
%   'spmversion'  = string, version of SPM to be used (default = 'spm12')
%
% The supported dataformats are
%   'afni_head'/'afni_brik'      uses AFNI code
%   'analyze_img'/'analyze_hdr'  uses SPM code
%   'analyze_old'                uses Darren Webber's code
%   'asa_mri'
%   'ctf_mri'
%   'ctf_mri4'
%   'ctf_svl'
%   'dicom'                      uses FreeSurfer code
%   'dicom_old'                  uses MATLAB image processing toolbox code
%   'freesurfer_mgh'             uses FreeSurfer code
%   'freesurfer_mgz'             uses FreeSurfer code
%   'jnifti_jnii'
%   'jnifti_bnii'
%   'matlab'                     assumes a MATLAB *.mat file containing a struct
%   'minc'                       uses SPM, this requires SPM5 or older
%   'mrtrix_mif'                 uses mrtrix code
%   'neuromag_fif'               uses MNE toolbox
%   'neuromag_fif_old'           uses meg-pd toolbox
%   'nifti'                      uses FreeSurfer code
%   'nifti_spm'                  uses SPM
%   'yokogawa_mri'
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
% If you have a series of DICOM files, please provide the name of any of the files in
% the series (e.g. the first one). The files corresponding to the whole volume will
% be found automatically.
%
% The output MRI may have a homogenous transformation matrix that converts the
% coordinates of each voxel (in xgrid/ygrid/zgrid) into head coordinates.
%
% If the input file is a 4D nifti, and you wish to load in just a subset of the
% volumes (e.g. due to memory constraints), you should use as dataformat 'nifti_spm',
% which uses the optional key-value pair 'volumes' = vector, with the indices of the
% to-be-read volumes, the order of the indices is ignored, and the volumes will be
% sorted according to the numeric indices, i.e. [1:10] yields the same as [10:-1:1]
%
% See also FT_DATATYPE_VOLUME, FT_WRITE_MRI, FT_READ_DATA, FT_READ_HEADER, FT_READ_EVENT

% Copyright (C) 2008-2022, Robert Oostenveld & Jan-Mathijs Schoffelen
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
spmversion  = ft_getopt(varargin, 'spmversion');

% use the version that is on the path, or default to spm12
if ~ft_hastoolbox('spm') && isempty(spmversion)
  spmversion = 'spm12';
end

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

if strcmp(dataformat, 'compressed') || (strcmp(dataformat, 'freesurfer_mgz') && ispc) || any(filetype_check_extension(filename, {'gz', 'zip', 'tar', 'tgz'}))
  % the file is compressed, unzip on the fly,
  % -freesurfer mgz files get special treatment only on a pc
  % -compressed AFNI BRIKS also need the HEAD copied over to the temp dir
  
  filename_old = filename;
  filename    = inflate_file(filename_old);
  if strcmp(dataformat, 'freesurfer_mgz')
    filename_old = filename;
    filename     = [filename '.mgh'];
    movefile(filename_old, filename);
  end
  if ~strcmp(dataformat, 'nifti_spm')
    % replace it with the filetype's default format, but don't overwrite in
    % case dataformat was nifti_spm
    dataformat = ft_filetype(filename);
  end
  if strcmp(dataformat, 'afni_brik')
    [p, f, e] =fileparts(filename_old);
    copyfile(fullfile(p, strrep(f, 'BRIK', 'HEAD')), strrep(filename, 'BRIK', 'HEAD'));
  end
  inflated = true;
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
      ft_hastoolbox(spmversion, 1);
    end
    volumes = ft_getopt(varargin, 'volumes', []);

    % use the functions from SPM
    hdr = spm_vol_nifti(filename);
    if isempty(volumes)
      img = double(hdr.private.dat);
    else
      volumes = sort(intersect(volumes, 1:size(hdr.private.dat,4)));
      img = double(hdr.private.dat(:,:,:,volumes));
    end
    transform = hdr.mat;
    unit = 'mm';

    try
      % The nifti header allows three methods for specifying the coordinates of the
      % voxels. For two of them (sform and qform), the header contains a numerical
      % code from 0 to 4 that specifies the coordinate system. SPM8 and SPM12 use a
      % field in the private header of the nifti object to code this.
      switch hdr.private.mat_intent
        case 'UNKNOWN'
          coordsys = 'unknown';
        case 'Scanner'
          coordsys = 'scanras';
        case 'Aligned'
          coordsys = 'aligned';
        case 'Talairach'
          coordsys = 'tal';
        case 'MNI152'
          coordsys = 'mni152';
      end

      % We cannot trust it yet, see https://github.com/fieldtrip/fieldtrip/issues/1879
      ft_notice('the coordinate system appears to be ''%s''\n', coordsys);
      clear coordsys
    end


  case {'analyze_img' 'analyze_hdr'}
    if ~(hasspm8 || hasspm12)
      fprintf('the SPM8 or newer toolbox is required to read analyze files\n');
      ft_hastoolbox(spmversion, 1);
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
    img = flip(img, 1);

    transform      = diag(hdr.dime.pixdim(2:4));
    transform(4,4) = 1;

  case {'afni_brik' 'afni_head'}
    % needs afni
    ft_hastoolbox('afni', 1);    % see http://afni.nimh.nih.gov/

    [err, hdr] = BrikInfo(filename);
    
    % check the precision of the data, and if scaling is required. If the precision is other than float, 
    %and no scaling is required, then return the data in its native precision, let the low level code take
    % care of that
    if any(hdr.BRICK_FLOAT_FACS~=0)
      opts.OutPrecision = '';
    else
      opts.OutPrecision = '*';
      opts.Scale = 0;
    end
    [err, img, hdr, ErrMessage] = BrikLoad(filename, opts);
    if err
      ft_error('could not read AFNI file');
    end

    if isfield(hdr, 'ORIENT_SPECIFIC')
      [err, orient, flipvec] = AFNI_OrientCode(hdr.ORIENT_SPECIFIC);
      % FIXME, I don't understand why the orient vector needs to be like
      % this: it seems the opposite of what is reflected in the coordsys
      % (see below), but it seems to yield internally consistent results
    else
      % afni volume info
      orient = 'LPI'; % hope for the best
    end
        
    % origin and basis vectors in world space
    [unused, ix] = AFNI_Index2XYZcontinuous([0 0 0; eye(3)], hdr, orient);
    
    % basis vectors in voxel space
    e1 = ix(2,:) - ix(1,:);
    e2 = ix(3,:) - ix(1,:);
    e3 = ix(4,:) - ix(1,:);

    % change from base0 (afni) to base1 (SPM/Matlab)
    o = ix(1,:) - (e1+e2+e3);

    % create matrix
    transform = [e1;e2;e3;o]';
    transform = cat(1, transform, [0 0 0 1]);
    
    coordsys = lower(hdr.Orientation(:,2)');
    if contains(filename, 'TTatlas') || (isfield(hdr, 'TEMPLATE_SPACE') && ~isempty(hdr.TEMPLATE_SPACE))
      if isfield(hdr, 'TEMPLATE_SPACE') && ~isempty(hdr.TEMPLATE_SPACE)
        space = hdr.TEMPLATE_SPACE;
      else
        space = 'tal'; % accommodate the case when this is not specified in the hdr, make assumption
      end
      if startsWith(space, 'tt_') || startsWith(space, 'TT_')
        space = 'tal';
      elseif startsWith(space, 'mni')
        space = 'mni';
      elseif startsWith(space, 'tlrc')
        % according to the documentation tlrc is rather generic as a
        % specification of the space, but originally it meant tal.
        ft_warning('space ''tlrc'' might be ambiguous, here assuming the coordsys to be ''tal''');
        space = 'tal';
      end
      if ismember(space, {'tal' 'mni'})
        % xyz orientation should be RAS
        if ~strcmp(coordsys, 'ras')
          ft_warning('the template space suggests that the image is in %s coordinates, but the xyz orientation %s does not match this', space, coordsys);
          xxx2ras = true;
        else
          coordsys = space;
        end
      end
    end

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
    % coordinate system 2002 - Surface RAS coordinates (is mainly vertical shift, no rotation to 2001)
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
      % centers the coordinate system and switches from mm to m
      if (hdr.voxel_trans.from == 2001) && (hdr.voxel_trans.to == 5)
        % matlab_shift compensates for the different index conventions between C and MATLAB

        % the lines below is old code (prior to Jan 3, 2013) and only works with 1 mm resolution MRIs
        %   matlab_shift = [ 0 0 0 0.001; 0 0 0 -0.001; 0 0 0 0.001; 0 0 0 0];
        % transform transforms from 2001 to 5 and further to 4
        %   transform = transform\(hdr.voxel_trans.trans+matlab_shift);

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
    % this returns a right-handed volume with the transformation matrix stored in the file headers
    % see https://github.com/fieldtrip/website/pull/444

    % needs the freesurfer toolbox
    ft_hastoolbox('freesurfer', 1);
    [dcmdir, junk1, junk2] = fileparts(filename);
    if isempty(dcmdir)
      dcmdir = '.';
    end
    [img, transform,hdr, mr_params] = load_dicom_series(dcmdir,dcmdir,filename);
    transform = vox2ras_0to1(transform);
    coordsys  = 'scanras';
    unit      = 'mm';

  case 'dicom_old'
    % this uses the Image processing toolbox
    % the DICOM files represent a stack of slices, and possibly even multiple volumes
    orig = dicominfo(filename);
    dim(1) = orig.Rows;
    dim(2) = orig.Columns;

    % this works for the Siemens scanners at the FCDC
    [p, f] = fileparts(filename);
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
      ft_info('reading image data from ''%s''\n', filename);
      img(:,:,i) = dicomread(hdr(i));
    end

    % reorder and concatenate the slices
    [z, indx]   = sort(cell2mat({hdr.SliceLocation}));
    hdr = hdr(indx);
    img = img(:,:,indx);

    % construct a homgeneous transformation matrix that performs the scaling from voxels to mm
    transform = dicom2transform(hdr);
    coordsys  = 'dicom'; % identical to scanlps, see https://www.fieldtriptoolbox.org/faq/coordsys/#details-of-the-dicom-coordinate-system
    unit      = 'mm';

    % this makes the mapping of voxels to patient coordinates consistent with Horos
    img = permute(img, [2, 1, 3]);

  case {'nifti', 'nifti_gz', 'freesurfer_mgz', 'freesurfer_mgh'}
    ft_hastoolbox('freesurfer', 1);
    tmp = MRIread(filename);
    ndims = numel(size(tmp.vol));
    if ndims==3
      img = permute(tmp.vol, [2 1 3]);
      % FIXME although this is probably correct
      % see the help of MRIread, anecdotally columns and rows seem to need a swap
      % in order to match the transform matrix (alternatively a row switch of the
      % latter can be done)
    elseif ndims==4
      img = permute(tmp.vol, [2 1 3 4]);
    end
    hdr = rmfield(tmp, 'vol');
    transform = hdr.vox2ras1;
    unit = 'mm';

    if isfield(hdr, 'niftihdr')
      % The nifti header allows three methods for specifying the coordinates of the
      % voxels. For two of them (sform and qform), the header contains a numerical
      % code that specifies the coordinate system.
      coordsys = {'unknown', 'scanras', 'aligned', 'tal', 'mni152'}; % corresponding to 0, 1, 2, 3, 4
      if isequal(hdr.vox2ras0, hdr.niftihdr.sform)
        coordsys = coordsys{hdr.niftihdr.sform_code + 1};
      elseif isequal(hdr.vox2ras0, hdr.niftihdr.qform)
        coordsys = coordsys{hdr.niftihdr.qform_code + 1};
      else
        coordsys = 'unknown';
      end

      % We cannot trust it yet, see https://github.com/fieldtrip/fieldtrip/issues/1879
      ft_notice('the coordinate system appears to be ''%s''\n', coordsys);
      clear coordsys
    end

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
    isfixel = numel(tmp.dim==3) && tmp.dim(3)==1;

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

  case {'openjdata_jnii' 'openjdata_bnii'}
    % this depends on two external toolboxes
    ft_hastoolbox('jsonlab', 1);
    ft_hastoolbox('jnifti', 1);

    jnii = loadjnifti(filename);

    mri.hdr     = jnii.NIFTIHeader;
    mri.anatomy = jnii.NIFTIData;
    mri.dim     = jnii.NIFTIHeader.Dim;
    mri.unit    = jnii.NIFTIHeader.Unit.L; % units of length

    % see https://brainder.org/2012/09/23/the-nifti-file-format/
    coordsys = {'scanras', 'aligned', 'tal', 'mni152'};

    if jnii.NIFTIHeader.SForm>0
      mri.coordsys  = coordsys{jnii.NIFTIHeader.SForm};
      mri.transform = [
        jnii.NIFTIHeader.Affine
        0 0 0 1
        ];

    elseif jnii.NIFTIHeader.QForm>0
      mri.coordsys = coordsys{jnii.NIFTIHeader.QForm};

      % this is adapted from freesurfer/load_nifti_hdr.m
      b = jnii.NIFTIHeader.Quatern.b;
      c = jnii.NIFTIHeader.Quatern.c;
      d = jnii.NIFTIHeader.Quatern.d;
      x = jnii.NIFTIHeader.QuaternOffset.x;
      y = jnii.NIFTIHeader.QuaternOffset.y;
      z = jnii.NIFTIHeader.QuaternOffset.z;

      a = 1.0 - (b*b + c*c + d*d);
      if(abs(a) < 1.0e-7)
        a = 1.0 / sqrt(b*b + c*c + d*d);
        b = b*a;
        c = c*a;
        d = d*a;
        a = 0.0;
      else
        a = sqrt(a);
      end

      r11 = a*a + b*b - c*c - d*d;
      r12 = 2.0*b*c - 2.0*a*d;
      r13 = 2.0*b*d + 2.0*a*c;
      r21 = 2.0*b*c + 2.0*a*d;
      r22 = a*a + c*c - b*b - d*d;
      r23 = 2.0*c*d - 2.0*a*b;
      r31 = 2.0*b*d - 2*a*c;
      r32 = 2.0*c*d + 2*a*b;
      r33 = a*a + d*d - c*c - b*b;

      if(jnii.NIFTIHeader.VoxelSize(1) < 0.0)
        r13 = -r13;
        r23 = -r23;
        r33 = -r33;
      end

      R = [r11 r12 r13; r21 r22 r23; r31 r32 r33];
      S = diag(jnii.NIFTIHeader.VoxelSize(2:4));
      T = [x y z]';
      mri.transform = [R*S T; 0 0 0 1];

    else
      mri.coordsys = 'unknown';
      % Method 1
      mri.transform      = diag(jnii.NIFTIHeader.VoxelSize(2:4));
      mri.transform(4,4) = 1;
    end

    % ensure that this is double precision and not uint8
     mri.transform = double(mri.transform);

    % these are already part of the output structure and should not be reassigned
    clear coordsys transform unit

  otherwise
    ft_error('unrecognized filetype ''%s'' for ''%s''', dataformat, filename);
end

if exist('img', 'var')
  % determine the size of the volume in voxels
  mri.dim = size(img);
  if numel(mri.dim)>3
    mri.dim = mri.dim(1:3); % dim should be the spatial 3D dim
  end
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

if exist('transform', 'var')
  % store the homogeneous transformation matrix if present
  mri.transform = transform;
end

if exist('unit', 'var')
  % determine the geometrical units in which it is expressed
  mri.unit = unit;
else
  % estimate the units from the data
  mri = ft_determine_units(mri);
end

if exist('coordsys', 'var')
  % add a descriptive label for the coordinate system
  mri.coordsys = coordsys;
end

if exist('xxx2ras', 'var') && xxx2ras==true
  % this is needed for AFNI formatted data, where the created voxels-to-world
  % mapping matrix is diagonal for the 3x3 rotation part (i.e. ijk should
  % be ras, in order for the tal/mni coordsys to make sense
  mri = ft_convert_coordsys(mri, 'ras', 0);
  mri.coordsys = space;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bool = isrighthanded(orient)
bool = ismember(orient, {'ALS' 'RAS' 'PRS' 'LPS' 'SAL' 'SRA' 'SPR' 'SLP'});
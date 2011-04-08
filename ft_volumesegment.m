function [segment] = ft_volumesegment(cfg, mri)

% FT_VOLUMESEGMENT segments an anatomical MRI into gray matter, white matter,
% and cerebro-spinal fluid compartments. It can also be used to create a
% binary mask, either of the segmented volumes, or of the 'raw' anatomical
% data.
%
% This function uses the SPM8 toolbox, see http://www.fil.ion.ucl.ac.uk/spm/
%
% Use as
%   [segment] = ft_volumesegment(cfg, mri)
%
% The input arguments are a configuration structure (see below) and an
% anatomical MRI structure. Instead of an MRI structure, you can also
% specify a string with a filename of an MRI file. You can also provide an
% already segmented volume in the input for the purpose of creating a
% binary mask of the inner surface of the skull.
%
% The configuration options are
%   cfg.spmversion  = 'spm8' (default) or 'spm2'
%   cfg.template    = filename of the template anatomical MRI (default is the 'T1.nii' 
%                     (spm8) or 'T1.mnc' (spm2) in the (spm-directory)/templates/)
%   cfg.name        = string for output filename
%   cfg.write       = 'no' or 'yes' (default = 'no'),
%                     writes the segmented volumes to SPM compatible analyze (spm2),
%                     or nifti (spm8) files,
%                     with the suffix (spm2)
%                     _seg1, for the gray matter segmentation
%                     _seg2, for the white matter segmentation
%                     _seg3, for the csf segmentation
%                     or with the prefix (spm8)
%                     c1, for the gray matter segmentation
%                     c2, for the white matter segmentation
%                     c3, for the csf segmentation
%                   
%   cfg.smooth      = 'no' (default) or the FWHM of the gaussian kernel in voxels
%   cfg.threshold   = 'no' or a threshold value which is used to threshold
%                     the data in order to create a volumetric mask (see below). 
%   cfg.segment     = 'yes' (default) or 'no'
%
% Example use:
%
%   segment = ft_volumesegment([], mri) will segment the anatomy and will output
%               the segmentation result as 3 probabilistic masks in 
%               segment.gray/.white/.csf
%
%   cfg = [];
%   cfg.smooth    = 5;
%   cfg.threshold = 0.5;
%   segment = ft_volumesegment(cfg, mri) will segment the anatomy and will
%               output not only the probabilistic masks, but also a binary
%               segment.brainmask, specifying the inside of the skull
%
%
%   cfg = [];
%   cfg.smooth    = 5;
%   cfg.threshold = 0.1;
%   cfg.segment   = 'no';
%   mri = ft_volumesegment(cfg, mri) will not segment the anatomy but only
%               smoothes and thresholds the anatomy in order to create a
%               binary mask mri.scalpmask, specifying the surface of the
%               scalp
%
% For the segmentation to work, the coordinate frame of the input MRI has to
% be approximately aligned to the template. For this, a homogeneous
% transformation matrix is used, which makes the assumption that the
% template mri is defined in SPM/MNI-coordinates:
%   x-axis pointing to the right ear
%   y-axis along the acpc-line
%   z-axis pointing to the top of the head
%   origin in the anterior commissure.
% Note that the segmentation only works if the template MRI is in SPM
% coordinates.
% 
% If the input mri is a string pointing to a CTF *.mri file, the
% x-axis is assumed to point to the nose, and the origin is assumed
% to be on the interauricular line.
%
% If the input mri is a string pointing to another fileformat, or a
% structure containing an anatomical MRI in Matlab memory, the user will
% be asked about the axis-definition and the origin of the coordinate system.
%
% As a second step, the segmentation is performed, using the
% default parameters from SPM. The output volume is in the original
% coordinate-frame.
% 
% As a third and optional step, you can perform a smoothing of the segmented
% volumes.
%
% To facilitate data-handling and distributed computing with the peer-to-peer
% module, this function has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_READ_MRI

% undocumented options
%   cfg.keepintermediate = 'yes' or 'no'

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

ft_defaults

cfg = ft_checkconfig(cfg, 'trackconfig', 'on');
cfg = ft_checkconfig(cfg, 'deprecated',  'coordinates');

%% ft_checkdata see below!!! %%

% set the defaults
cfg.segment          = ft_getopt(cfg, 'segment',          'yes');
cfg.smooth           = ft_getopt(cfg, 'smooth',           'no');
cfg.spmversion       = ft_getopt(cfg, 'spmversion',       'spm8');
cfg.write            = ft_getopt(cfg, 'write',            'no');
cfg.threshold        = ft_getopt(cfg, 'threshold',        'no');
cfg.keepintermediate = ft_getopt(cfg, 'keepintermediate', 'no');
cfg.coordsys         = ft_getopt(cfg, 'coordsys',         '');
cfg.units            = ft_getopt(cfg, 'units',            '');
cfg.inputfile        = ft_getopt(cfg, 'inputfile',        []);
cfg.outputfile       = ft_getopt(cfg, 'outputfile',       []);

if ~isfield(cfg,'name') 
  if ~strcmp(cfg.write,'yes')
    tmp = tempname;
    cfg.name = tmp;
  else
    error('you must specify the output filename in cfg.name');
  end
end 

% check if the required spm is in your path:
if strcmpi(cfg.spmversion, 'spm2'),
  ft_hastoolbox('SPM2',1);
elseif strcmpi(cfg.spmversion, 'spm8'),
  ft_hastoolbox('SPM8',1);
end

% get the names of the templates for the segmentation
if ~isfield(cfg, 'template'),
  spmpath      = spm('dir');
  if strcmpi(cfg.spmversion, 'spm8'), cfg.template = [spmpath,filesep,'templates',filesep,'T1.nii']; end
  if strcmpi(cfg.spmversion, 'spm2'), cfg.template = [spmpath,filesep,'templates',filesep,'T1.mnc']; end
end

hasdata      = (nargin>1);
hasinputfile = ~isempty(cfg.inputfile);
if hasdata && hasinputfile
  error('cfg.inputfile should not be used in conjunction with giving input data to this function');
elseif hasinputfile
  % the input data should be read from file
  mri = loadvar(cfg.inputfile, 'mri');
elseif hasdata
  if ischar(mri),
    % read the anatomical MRI data from file
    filename = mri;
    fprintf('reading MRI from file\n');
    mri = ft_read_mri(filename);
    if ft_filetype(filename, 'ctf_mri') && isempty(cfg.coordsys)
      % based on the filetype assume that the coordinates correspond with CTF convention
      cfg.coordsys = 'ctf';
    end
  end
else
  error('neither a data structure, nor a cfg.inputfile is provided');
end

% check if the input data is valid for this function
mri = ft_checkdata(mri, 'datatype', 'volume', 'feedback', 'yes');

% check whether the input data contains already segmented data
hassegment = 0;
hasanatomy = 0;
hasbrain   = 0;
if isfield(mri, 'anatomy'),
  hasanatomy = 1;
end

if isfield(mri, 'gray') && isfield(mri, 'white') && isfield(mri, 'csf'),
  hassegment = 1;
  if strcmp(cfg.segment, 'yes')
    fprintf('the input data already contains a segmentation of gray/white/csf matter, no segmentation will be performed\n');
    cfg.segment = 'no';
  end
end

if strcmp(cfg.segment, 'yes')
  
  % ensure that the data has interpretable units and that the coordinate
  % system is in approximate spm space
  if ~isfield(mri, 'unit'),     mri.unit     = cfg.units;    end
  if ~isfield(mri, 'coordsys'), mri.coordsys = cfg.coordsys; end
  
  % remember the original transformation matrix coordinate system
  original.transform = mri.transform;
  original.coordsys  = mri.coordsys; 
  
  mri = ft_convert_units(mri,    'mm');
  mri = ft_convert_coordsys(mri, 'spm');
  
  % flip and permute the 3D volume itself, so that the voxel and 
  % headcoordinates approximately correspond this improves the convergence
  % of the segmentation algorithm
  [mri,permutevec,flipflags] = align_ijk2xyz(mri);
  
  Va = ft_write_volume([cfg.name,'.img'], mri.anatomy, 'transform', mri.transform, 'spmversion', cfg.spmversion);

  % spm is quite noisy, prevent the warnings from displaying on screen
  % warning off;

  if strcmpi(cfg.spmversion, 'spm2'),
    % set the spm segmentation defaults (from /opt/spm2/spm_defaults.m script)
    defaults.segment.estimate.priors = str2mat(...
      fullfile(spm('Dir'),'apriori','gray.mnc'),...
      fullfile(spm('Dir'),'apriori','white.mnc'),...
      fullfile(spm('Dir'),'apriori','csf.mnc'));
    defaults.segment.estimate.reg    = 0.01;
    defaults.segment.estimate.cutoff = 30;
    defaults.segment.estimate.samp   = 3;
    defaults.segment.estimate.bb     =  [[-88 88]' [-122 86]' [-60 95]'];
    defaults.segment.estimate.affreg.smosrc = 8;
    defaults.segment.estimate.affreg.regtype = 'mni';
    %defaults.segment.estimate.affreg.weight = fullfile(spm('Dir'),'apriori','brainmask.mnc'); 
    defaults.segment.estimate.affreg.weight = '';
    defaults.segment.write.cleanup   = 1;
    defaults.segment.write.wrt_cor   = 1;
    
    flags = defaults.segment;

    % perform the segmentation
    fprintf('performing the segmentation on the specified volume\n');
    spm_segment(Va,cfg.template,flags);
    Vtmp = spm_vol({[cfg.name,'_seg1.img'];...
                    [cfg.name,'_seg2.img'];...
                    [cfg.name,'_seg3.img']});

    % read the resulting volumes
    for j = 1:3
      vol = spm_read_vols(Vtmp{j});
      Vtmp{j}.dat = vol;
      V(j) = struct(Vtmp{j});
    end

    % keep or remove the files according to the configuration
    if strcmp(cfg.keepintermediate,'no'),
      delete([cfg.name,'.img']);
      delete([cfg.name,'.hdr']);
      delete([cfg.name,'.mat']);
    end
    if strcmp(cfg.write,'no'),
       delete([cfg.name,'_seg1.hdr']);
       delete([cfg.name,'_seg2.hdr']);
       delete([cfg.name,'_seg3.hdr']);
       delete([cfg.name,'_seg1.img']);
       delete([cfg.name,'_seg2.img']);
       delete([cfg.name,'_seg3.img']);
       delete([cfg.name,'_seg1.mat']);
       delete([cfg.name,'_seg2.mat']);
       delete([cfg.name,'_seg3.mat']);
    elseif strcmp(cfg.write,'yes'),
      for j = 1:3
        % put the original transformation-matrix in the headers
        V(j).mat = original.transform;
        % write the updated header information back to file ???????
        V(j) = spm_create_vol(V(j));
      end
    end

  elseif strcmpi(cfg.spmversion, 'spm8'),
    
    fprintf('performing the segmentation on the specified volume\n');
    if isfield(cfg, 'tpm')
      px.tpm   = cfg.tpm;
      p        = spm_preproc(Va, px);
    else
      p        = spm_preproc(Va);
    end
    [po,pin] = spm_prep2sn(p);
    
    % I took these settings from a batch
    opts     = [];
    opts.GM  = [0 0 1];
    opts.WM  = [0 0 1];
    opts.CSF = [0 0 1];
    opts.biascor = 1;
    opts.cleanup = 0;
    spm_preproc_write(po, opts);
     
    [pathstr,name,ext] = fileparts(cfg.name);
    Vtmp = spm_vol({fullfile(pathstr,['c1',name,'.img']);...
                    fullfile(pathstr,['c2',name,'.img']);...
                    fullfile(pathstr,['c3',name,'.img'])});

    % read the resulting volumes
    for j = 1:3
      vol = spm_read_vols(Vtmp{j});
      Vtmp{j}.dat = vol;
      V(j) = struct(Vtmp{j});
    end

    % keep or remove the files according to the configuration
    if strcmp(cfg.keepintermediate,'no'),
      delete([cfg.name,'.img']);
      delete([cfg.name,'.hdr']);
      if exist([cfg.name,'.mat'], 'file'), 
        delete([cfg.name,'.mat']);
      end %does not always exist
    end
    
    % keep the files written to disk or remove them
    % FIXME check whether this works at all
    if strcmp(cfg.write,'no'),
       delete(fullfile(pathstr,['c1',name,'.hdr'])); %FIXME this may not be needed in spm8
       delete(fullfile(pathstr,['c1',name,'.img']));
       delete(fullfile(pathstr,['c2',name,'.hdr']));
       delete(fullfile(pathstr,['c2',name,'.img']));
       delete(fullfile(pathstr,['c3',name,'.hdr']));
       delete(fullfile(pathstr,['c3',name,'.img']));
       delete(fullfile(pathstr,['m',name,'.hdr']));
       delete(fullfile(pathstr,['m',name,'.img']));
    elseif strcmp(cfg.write,'yes'),
      for j = 1:3
        % put the original transformation-matrix in the headers
        V(j).mat = original.transform;
        % write the updated header information back to file ???????
        V(j) = spm_create_vol(V(j));
      end
    end
    
  end

  % collect the results
  segment.dim       = size(V(1).dat);
  
  segment.dim       = segment.dim(:)';    % enforce a row vector
  segment.transform = original.transform; % use the original transform
  segment.coordsys  = original.coordsys;  % use the original coordsys
  if isfield(mri, 'unit')
    segment.unit = mri.unit;
  end
  segment.gray      = V(1).dat;
  if length(V)>1, segment.white     = V(2).dat; end
  if length(V)>2, segment.csf       = V(3).dat; end

  % flip the volumes back according to the changes introduced by align_ijk2xyz
  for k = 1:3
    if flipflags(k)
      segment.gray = flipdim(segment.gray, k);
      if isfield(segment, 'white'), segment.white = flipdim(segment.white, k); end
      if isfield(segment, 'csf'),   segment.csf   = flipdim(segment.csf, k);   end
    end
  end

  if ~all(permutevec == [1 2 3])
    segment.gray = ipermute(segment.gray, permutevec);
    if isfield(segment, 'white'), segment.white = ipermute(segment.white, permutevec); end
    if isfield(segment, 'csf'),   segment.csf   = ipermute(segment.csf,   permutevec); end
    segment.dim  = size(segment.gray);
  end

  hassegment = 1;
else
  % rename the data
  segment = mri;
  clear mri;  
end

% prepare for optional thresholding
if ~strcmp(cfg.threshold, 'no'),
  if hassegment
    % combine the single image segmentation from the three probabilistic
    % tissue segmentations for csf, white and gray matter
    segment.brain = segment.gray + segment.white + segment.csf;
    hasbrain = 1;
  end
end

% do optional smoothing
if ~strcmp(cfg.smooth, 'no'),
  fprintf('smoothing with a %d-voxel FWHM kernel\n',cfg.smooth);

  if hassegment
    spm_smooth(segment.gray,  segment.gray,  cfg.smooth); % smooth the gray matter
    spm_smooth(segment.white, segment.white, cfg.smooth); % smooth the white matter
    spm_smooth(segment.csf,   segment.csf,   cfg.smooth); % smooth the csf compartment
  end
  if hasanatomy & strcmp(cfg.segment,'no')
    spm_smooth(segment.anatomy, segment.anatomy, cfg.smooth); 
  else
    segment.anatomy = mri.anatomy;
    spm_smooth(segment.anatomy, segment.anatomy, cfg.smooth); 
  end
  if hasbrain
    spm_smooth(segment.brain, segment.brain, cfg.smooth);
  end
end

% do optional thresholding
if ~strcmp(cfg.threshold, 'no'),
  fprintf('thresholding the data at a relative threshold of %0.3f\n',cfg.threshold);
  
  if hassegment
  end
  if hasanatomy
    % create scalp mask by taking the negative of the brain, thus ensuring
    % that no holes are within the 'head compartment' and do a two-pass 
    % approach to eliminate potential vitamin E capsules etc.
    segment.scalpmask = double(segment.anatomy>(cfg.threshold*max(segment.anatomy(:))));
    [tmp, N]          = spm_bwlabel(segment.scalpmask, 6);
    for k = 1:N
      n(k,1) = sum(tmp(:)==k);
    end
    segment.scalpmask = double(tmp~=find(n==max(n))); clear tmp;
    [tmp, N]          = spm_bwlabel(segment.scalpmask, 6);
    for k = 1:N
      m(k,1) = sum(tmp(:)==k);
    end
    segment.scalpmask = tmp~=find(m==max(m)); clear tmp;
  end
  if hasbrain
    % create brain mask
    segment.brainmask = segment.brain>(cfg.threshold*max(segment.brain(:)));
    segment = rmfield(segment, 'brain');
  end
    
end

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add version information to the configuration
cfg.version.name = mfilename('fullpath');
cfg.version.id = '$Id$';

% add information about the Matlab version used to the configuration
cfg.version.matlab = version();

% remember the configuration details of the input data
if isfield(segment, 'cfg'),
  cfg.previous = segment.cfg;
end

% remember the exact configuration details in the output 
segment.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'segment', segment); % use the variable name "segment" in the output file
end


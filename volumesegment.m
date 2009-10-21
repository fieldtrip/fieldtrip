function [segment] = volumesegment(cfg, mri)

% VOLUMESEGMENT segments an anatomical MRI into gray matter, white matter,
% and cerebro-spinal fluid compartments.
%
% This function the SPM2 toolbox, see http://www.fil.ion.ucl.ac.uk/spm/
%
% Use as
%   [segment] = volumesegment(cfg, mri)
%
% The input arguments are a configuration structure (see below) and an
% anatomical MRI structure. Instead of an MRI structure, you can also
% specify a string with a filename of an MRI file.
%
% The configuration options are
%   cfg.template    = filename of the template anatomical MRI (default is '/home/common/matlab/spm2/templates/T1.mnc')
%   cfg.name        = string for output filename
%   cfg.write       = 'no' or 'yes' (default = 'no'),
%                     writes the segmented volumes to SPM2 compatible analyze-file with the suffix
%                     _seg1, for the gray matter segmentation
%                     _seg2, for the white matter segmentation
%                     _seg3, for the csf segmentation
%   cfg.smooth      = 'no' or the FWHM of the gaussian kernel in voxels (default = 'no')
%   cfg.coordinates = 'spm, 'ctf' or empty for interactive (default = [])
%
% As the first step the coordinate frame of the input MRI has to
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

% undocumented options
%   cfg.keepintermediate = 'yes' or 'no'
%   cfg.segment          = 'yes' or 'no'

% $Log: volumesegment.m,v $
% Revision 1.14  2009/07/14 07:27:30  roboos
% replaced read_fcdc_mri with read_mri to avoid warning
%
% Revision 1.13  2009/01/20 13:01:31  sashae
% changed configtracking such that it is only enabled when BOTH explicitly allowed at start
% of the fieldtrip function AND requested by the user
% in all other cases configtracking is disabled
%
% Revision 1.12  2008/11/21 13:56:12  sashae
% added call to checkconfig at start and end of function
%
% Revision 1.11  2008/09/22 20:17:44  roboos
% added call to fieldtripdefs to the begin of the function
%
% Revision 1.10  2008/09/17 14:53:35  roboos
% removed fixvolume (and underlying grid2transform), not needed any more because checkdata has the possibility of converting a pos to a transform
%
% Revision 1.9  2007/07/31 13:00:05  jansch
% minor change in save-directory for temporary files
%
% Revision 1.8  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.7  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.6  2006/08/01 10:06:35  marsie
% fixed bug in using spm_smooth
%
% Revision 1.5  2006/07/27 08:28:46  roboos
% use spm_smooth instead of spm_conv, updated documentation
%
% Revision 1.4  2006/03/06 09:31:43  roboos
% changed the default template from /opt to /home/common
%
% Revision 1.3  2006/02/24 16:47:19  roboos
% use fixvolume instead of grid2transform
%
% Revision 1.2  2006/01/30 13:47:31  roboos
% change in whitespace
%
% Revision 1.1  2006/01/05 12:58:20  roboos
% This function (VOLUMExxx) replaces a function with the name xxxVOLUME.
% The fields xgrid/ygrid/zgrid are removed (this is from now on handled by
% grid2transform and the VOLUMExxx function should only work on volumes that
% are described using a transformation matrix).
% Writing of spm/analyze volumes is handled by private/volumewrite_spm.
%
% Revision 1.10  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.9  2005/05/17 17:50:38  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.8  2005/04/25 16:06:23  roboos
% added the documentation for cfg.smooth and cfg.coordinates to the help
% fixed small bugsmooth due to which it would never smooth
%
% Revision 1.7  2005/03/03 10:48:16  roboos
% Cleaned up the function in total. It now uses a different (hopefully
% more robust) initial pre-alignment of the input volume with the SPM
% template. Also the type-detection (i.e. either SPM or CTF) of the input
% volume has been changed.
%
% Revision 1.6  2004/10/14 07:01:33  jansch
% included preprocessing of input-mri, to get an approximate aligmnent with
% the spm-templates used, prior to the segmentation.
%
% Revision 1.5  2004/10/13 20:28:44  jansch
% added the computation of a transformation-matrix, in order to align the 
% to-be-segmented volume  with the template and the segmentation-masks
%
% Revision 1.4  2004/09/22 10:46:31  roboos
% ensure that the coordinate transformation is from 1:N to headcoordinates
%
% Revision 1.3  2004/09/20 13:29:42  roboos
% fixed some small bugs, changed layout of help and comments
%
% Revision 1.2  2004/09/17 13:47:02  jansch
% inserted some changes suggested by robert
%
% Revision 1.1  2004/09/17 11:23:59  jansch
% first implementation
%

fieldtripdefs

cfg = checkconfig(cfg, 'trackconfig', 'on');

%% checkdata see below!!! %%

% check if spm2 is in your path:
hasspm = (exist('spm_vol') & exist('spm_write_vol') & exist('spm_segment'));
if ~hasspm
  error('the SPM2 toolbox is not installed, see: http://www.fil.ion.ucl.ac.uk/spm/');
end

% set the defaults
if ~isfield(cfg,'segment'),         cfg.segment = 'yes';                                        end;
if ~isfield(cfg,'smooth'),          cfg.smooth = 'no';                                          end;
if ~isfield(cfg,'template'),        cfg.template = '/home/common/matlab/spm2/templates/T1.mnc'; end;
if ~isfield(cfg,'write'),           cfg.write = 'no';                                           end;
if ~isfield(cfg,'keepintermediate'),cfg.keepintermediate = 'no';                                end;
if ~isfield(cfg,'coordinates'),     cfg.coordinates = [];                                       end;

if ~isfield(cfg,'name') 
  if ~strcmp(cfg.write,'yes')
    tmp = tempname;
    %[path,file] = fileparts(tmp);
    %cfg.name = file;
    cfg.name = tmp;
  else
    error('you must specify the output filename in cfg.name');
  end
end 

if ischar(mri),
  % read the anatomical MRI data from file   
  filename = mri;
  fprintf('reading MRI from file\n');
  mri = read_mri(filename);
  if filetype(filename, 'ctf_mri') && isempty(cfg.coordinates)
    % based on the filetype assume that the coordinates correspond with CTF convention
    cfg.coordinates = 'ctf';
  end
end

% check if the input data is valid for this function
mri = checkdata(mri, 'datatype', 'volume', 'feedback', 'yes');

% remember the original transformation matrix
original.transform = mri.transform;
original.nosedir   = [];
original.origin    = [];

if isempty(cfg.coordinates)
  % determine in which direction the nose is pointing
  while isempty(original.nosedir)
    response = input('Is the nose pointing in the positive X- or Y-direction? [x/y] ','s');
    if strcmp(response, 'x'),
      original.nosedir = 'positivex';
    elseif strcmp(response, 'y'),
      original.nosedir = 'positivey';
    end
  end

  % determine where the origin is
  while isempty(original.origin)
    response = input('Is the origin on the Anterior commissure or between the Ears? [a/e] ','s');
    if strcmp(response, 'e'),
      original.origin = 'interauricular';
    elseif strcmp(response, 'a'),
      original.origin = 'ant_comm';
    end
  end

  % determine the coordinatesystem of the input MRI volume
  if strcmp(original.origin, 'interauricular') && strcmp(original.nosedir, 'positivex')
    cfg.coordinates = 'ctf';
  elseif strcmp(original.origin, 'ant_comm') && strcmp(original.nosedir, 'positivey')
    cfg.coordinates = 'spm';
  end
end

% ensure that the input MRI volume is approximately aligned with the SPM template
if strcmp(cfg.coordinates, 'ctf')
  fprintf('assuming CTF coordinates for input, i.e. positive X-axis towards nasion and Y-axis through ears\n');
  % flip, rotate and translate the CTF mri so that it approximately corresponds with SPM coordinates
  % this only involves a manipulation of the coordinate tarnsformation matrix
  mri = align_ctf2spm(mri);
  % also flip and permute the 3D volume itself, so that the voxel and headcoordinates approximately correspond
  % this seems to improve the convergence of the segmentation algorithm
  mri = align_ijk2xyz(mri);
elseif strcmp(cfg.coordinates, 'spm')
  fprintf('assuming that the input MRI is already approximately alligned with SPM coordinates\n');
  % nothing needs to be done
else
  error('cannot determine the (approximate) alignmenmt of the input MRI with the SPM template');
end

if strcmp(cfg.segment, 'yes')
  % convert and write the volume to an analyze format, so that it can be handled by spm
  Va = volumewrite_spm([cfg.name,'.img'], mri.anatomy, mri.transform);

  % spm is quite noisy, prevent the warnings from displaying on screen
  % warning off;

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

else
  % the volume is already segmented, put it in an SPM structure
  V     = mri;
  V.dat = V.anatomy;
  V.mat = V.transform;
  V     = rmfield(V,'anatomy');
end

if ~strcmp(cfg.smooth, 'no'),
  fprintf('smoothing with a %d-pixel FWHM kernel\n',cfg.smooth);
                  spm_smooth(V(1).dat, V(1).dat, cfg.smooth);     % smooth the gray matter
  if length(V)>1, spm_smooth(V(2).dat, V(2).dat, cfg.smooth); end % smooth the white matter
  if length(V)>2, spm_smooth(V(3).dat, V(3).dat, cfg.smooth); end % smooth the csf
end

% collect the results
segment.dim       = size(V(1).dat);
segment.transform = original.transform;           % use the original transformation-matrix
segment.gray      = V(1).dat;
if length(V)>1, segment.white     = V(2).dat; end
if length(V)>2, segment.csf       = V(3).dat; end

% get the output cfg
cfg = checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i).name;
end
cfg.version.id = '$Id: volumesegment.m,v 1.14 2009/07/14 07:27:30 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = mri.cfg; end
% remember the exact configuration details in the output 
segment.cfg = cfg;


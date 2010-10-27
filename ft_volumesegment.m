function [segment] = ft_volumesegment(cfg, mri)

% FT_VOLUMESEGMENT segments an anatomical MRI into gray matter, white matter,
% and cerebro-spinal fluid compartments.
%
% This function uses the SPM2 toolbox, see http://www.fil.ion.ucl.ac.uk/spm/
%
% Use as
%   [segment] = ft_volumesegment(cfg, mri)
%
% The input arguments are a configuration structure (see below) and an
% anatomical MRI structure. Instead of an MRI structure, you can also
% specify a string with a filename of an MRI file.
%
% The configuration options are
%   cfg.spmversion  = 'spm8' (default) or 'spm2'
%   cfg.template    = filename of the template anatomical MRI (default is the 'T1.mnc' (spm2) or 'T1.nii' (spm8)
%                     in the (spm-directory)/templates/)
%   cfg.name        = string for output filename
%   cfg.write       = 'no' or 'yes' (default = 'no'),
%                     writes the segmented volumes to SPM compatible analyze-files,
%                     with the suffix (spm2)
%                     _seg1, for the gray matter segmentation
%                     _seg2, for the white matter segmentation
%                     _seg3, for the csf segmentation
%                     or with the prefix (spm8)
%                     c1, for the gray matter segmentation
%                     c2, for the white matter segmentation
%                     c3, for the csf segmentation
%                   
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
%
% See also FT_READ_MRI

% undocumented options
%   cfg.keepintermediate = 'yes' or 'no'
%   cfg.segment          = 'yes' or 'no'
%   cfg.inputfile        = one can specifiy preanalysed saved data as input
%   cfg.outputfile       = one can specify output as file to save to disk  

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

fieldtripdefs

cfg = ft_checkconfig(cfg, 'trackconfig', 'on');

%% ft_checkdata see below!!! %%


% set the defaults
if ~isfield(cfg,'segment'),         cfg.segment     = 'yes';                                     end;
if ~isfield(cfg,'smooth'),          cfg.smooth      = 'no';                                      end;
if ~isfield(cfg,'spmversion'),      cfg.spmversion  = 'spm8';                                    end;
if ~isfield(cfg,'write'),           cfg.write       = 'no';                                      end;
if ~isfield(cfg,'keepintermediate'),cfg.keepintermediate = 'no';                                 end;
if ~isfield(cfg,'coordinates'),     cfg.coordinates = [];                                        end;
if ~isfield(cfg, 'inputfile'),      cfg.inputfile   = [];   end
if ~isfield(cfg, 'outputfile'),     cfg.outputfile  = [];   end

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

% check if the required spm is in your path:
if strcmpi(cfg.spmversion, 'spm2'),
  ft_hastoolbox('SPM2',1);
elseif strcmpi(cfg.spmversion, 'spm8'),
  ft_hastoolbox('SPM8',1);
end

if ~isfield(cfg, 'template'),
  spmpath      = spm('dir');
  if strcmpi(cfg.spmversion, 'spm8'), cfg.template = [spmpath,filesep,'templates',filesep,'T1.nii']; end
  if strcmpi(cfg.spmversion, 'spm2'), cfg.template = [spmpath,filesep,'templates',filesep,'T1.mnc']; end
end

hasdata = (nargin>1);
if ~isempty(cfg.inputfile)
  % the input data should be read from file
  if hasdata
    error('cfg.inputfile should not be used in conjunction with giving input data to this function');
  else
    mri = loadvar(cfg.inputfile, 'data');
  end
end

if hasdata
  if ischar(mri),
    % read the anatomical MRI data from file
    filename = mri;
    fprintf('reading MRI from file\n');
    mri = ft_read_mri(filename);
    if ft_filetype(filename, 'ctf_mri') && isempty(cfg.coordinates)
      % based on the filetype assume that the coordinates correspond with CTF convention
      cfg.coordinates = 'ctf';
    end
  end
end

% check if the input data is valid for this function
mri = ft_checkdata(mri, 'datatype', 'volume', 'feedback', 'yes');

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
  fprintf('assuming that the input MRI is already approximately aligned with SPM coordinates\n');
  % nothing needs to be done
else
  error('cannot determine the (approximate) alignmenmt of the input MRI with the SPM template');
end

if strcmp(cfg.segment, 'yes')
  % convert and write the volume to an analyze format, so that it can be handled by spm
  Va = volumewrite_spm([cfg.name,'.img'], mri.anatomy, mri.transform, cfg.spmversion);

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
    p        = spm_preproc(Va);
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
      delete([cfg.name,'.mat']);
    end
    if strcmp(cfg.write,'no'),
       delete(fullfile(pathstr,['c1',name,'.hdr']));
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

% accessing this field here is needed for the configuration tracking
% by accessing it once, it will not be removed from the output cfg
cfg.outputfile;

% get the output cfg
cfg = ft_checkconfig(cfg, 'trackconfig', 'off', 'checksize', 'yes'); 

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i).name;
end
cfg.version.id = '$Id$';
% remember the configuration details of the input data
try, cfg.previous = mri.cfg; end
% remember the exact configuration details in the output 
segment.cfg = cfg;

% the output data should be saved to a MATLAB file
if ~isempty(cfg.outputfile)
  savevar(cfg.outputfile, 'data', segment); % use the variable name "data" in the output file
end



function [segment] = ft_volumesegment(cfg, mri)

% FT_VOLUMESEGMENT segments an anatomical MRI. The behaviour depends on the output requested. It can
% return probabilistic tissue maps of gray/white/csf compartments, a skull-stripped anatomy, or 
% binary masks representing the brain surface, skull, or scalp surface.
%
% Use as
%   segment = ft_volumesegment(cfg, mri)
% where the input mri should be a single anatomical volume that was for example read with 
% FT_READ_MRI. For the purpose of creating binary masks of the brain or of the skull, you can also 
% provide either the anatomical volume or the already segmented volume (with the probabilistic 
% tissue maps) as input.
%
% The configuration structure can contain
%   cfg.output         = string or cell-array of strings, see below (default = 'tpm')
%   cfg.spmversion     = string, 'spm2' or 'spm8' (default = 'spm8')
%   cfg.template       = filename of the template anatomical MRI (default = '/spm2/templates/T1.mnc'
%                        or '/spm8/templates/T1.nii')
%   cfg.name           = string for output filename
%   cfg.write          = 'no' or 'yes' (default = 'no'),
%                       writes the probabilistic tissue maps to SPM compatible analyze (spm2),
%                       or nifti (spm8) files, 
%                       with the suffix (spm2)
%                       _seg1, for the gray matter segmentation
%                       _seg2, for the white matter segmentation
%                       _seg3, for the csf segmentation
%                       or with the prefix (spm8)
%                       c1, for the gray matter segmentation
%                       c2, for the white matter segmentation
%                       c3, for the csf segmentation
%   cfg.brainsmooth    = 'no', or scalar, the FWHM of the gaussian kernel in voxels, (default = 5)
%   cfg.scalpsmooth    = 'no', or scalar, the FWHM of the gaussian kernel in voxels, (default = 5)
%   cfg.brainthreshold = 'no', or scalar, relative threshold value which is used to threshold the 
%                       tpm in order to create a volumetric brainmask (see below), (default = 0.5)
%   cfg.scalpthreshold = 'no', or scalar, relative threshold value which is used to threshold the 
%                       anatomical data in order to create a volumetric scalpmask (see below), 
%                       (default = 0.1)
%   cfg.downsample     = integer, amount of downsampling before segmentation 
%                       (default = 1; i.e., no downsampling)
%   cfg.coordsys       = string, specifying the coordinate system in which the anatomical data is 
%                       defined. This will be used if the input mri does not contain a 
%                       coordsys-field, (default = '', which results in the user being forced to 
%                       evaluate the coordinate system)
%   cfg.units          = the physical units in which the output will be expressed. (default = 'mm')
%
% The desired segmentation output is specified with cfg.output as a string or cell-array of strings
% and can contain
%   'tpm'         - tissue probability map for csf, white and gray matter
%   'brain'       - binary representation of the brain (including csf, white and gray matter) 
%   'skull'       - binary representation of the skull
%   'scalp'       - binary representation of the scalp
%   'skullstrip'  - anatomy with only the brain
%
% Example use:
%   cfg        = [];
%   segment    = ft_volumesegment(cfg, mri) will segment the anatomy and will output the 
%                segmentation result as 3 probabilistic masks in segment.gray, white and csf.
%
%   cfg.output = 'skullstrip';
%   segment    = ft_volumesegment(cfg, mri) will generate a skull-stripped anatomy based on a 
%                brainmask generated from the probabilistic tissue maps. The skull-stripped anatomy
%                is stored in the field segment.anatomy.
%
%
%   cfg.output = {'brain' 'scalp' 'skull'};
%   segment    = ft_volumesegment(cfg, mri) will produce a volume with 3 binary masks, representing
%                the brain surface, scalp surface, and skull.
%
%   cfg.output = {'scalp'};
%   segment    = ft_volumesegment(cfg, mri) will produce a volume with a binary mask (based on the 
%                anatomy), representing the border of the scalp surface (i.e., everything inside the 
%                surface is also included). Such representation of the scalp is produced faster, 
%                because it doesn't require to create the tissue probabilty maps before creating 
%                the mask. 
%
% For the SPM-based segmentation to work, the coordinate frame of the input MRI needs to be 
% approximately coregistered to the templates of the probabilistic tissue maps. The templates are 
% defined in SPM/MNI-space. FieldTrip attempts to do an automatic alignment based on the
% coordsys-field in the mri, and if this is not present, based on the coordsys-field in the cfg. If 
% none of them is specified the FT_DETERMINE_COORDSYS function is used to interactively assess the
% coordinate system in which the MRI is expressed.
%
% The template mri is defined in SPM/MNI-coordinates:
%   x-axis pointing to the right ear
%   y-axis along the acpc-line
%   z-axis pointing to the top of the head
%   origin in the anterior commissure.
% Note that the segmentation only works if the template MRI is in SPM coordinates.
% 
% If the input mri is a string pointing to a CTF *.mri file, the x-axis is assumed to point to the 
% nose, and the origin is assumed to be on the interauricular line. In this specific case, when 
% ft_read_mri is used to read in the mri, the coordsys field is automatically attached.
%
% To facilitate data-handling and distributed computing with the peer-to-peer module, this function 
% has the following options:
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat file on disk and/or
% the output data will be written to a *.mat file. These mat files should contain only a single 
% variable, corresponding with the input/output structure.
%
% See also FT_READ_MRI, FT_DETERMINE_COORDSYS, FT_PREPARE_HEADMODEL

% undocumented options
%   cfg.keepintermediate = 'yes' or 'no'

% Copyright (C) 2007-2012, Jan-Mathijs Schoffelen, Robert Oostenveld
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

revision = '$Id$';

% do the general setup of the function
ft_defaults
ft_preamble help
ft_preamble callinfo
ft_preamble trackconfig
ft_preamble loadvar mri

% this is not supported any more as of 26/10/2011
if ischar(mri),
  error('please use cfg.inputfile instead of specifying the input variable as a string');
end

% check if the input data is valid for this function
mri = ft_checkdata(mri, 'datatype', 'volume', 'feedback', 'yes');

% check if the input cfg is valid for this function
cfg = ft_checkconfig(cfg, 'renamed',  {'coordinates', 'coordsys'});

% set the defaults
cfg.output           = ft_getopt(cfg, 'output',           'tpm');
cfg.downsample       = ft_getopt(cfg, 'downsample',       1);
cfg.spmversion       = ft_getopt(cfg, 'spmversion',       'spm8');
cfg.write            = ft_getopt(cfg, 'write',            'no');
cfg.keepintermediate = ft_getopt(cfg, 'keepintermediate', 'no');
cfg.coordsys         = ft_getopt(cfg, 'coordsys',         '');
cfg.units            = ft_getopt(cfg, 'units',            '');
cfg.inputfile        = ft_getopt(cfg, 'inputfile',        []);
cfg.outputfile       = ft_getopt(cfg, 'outputfile',       []);

% set default for smooth and threshold
cfg.brainsmooth      = ft_getopt(cfg, 'brainsmooth',      ''); % see also below
cfg.scalpsmooth      = ft_getopt(cfg, 'scalpsmooth',      ''); % see also below
cfg.brainthreshold   = ft_getopt(cfg, 'brainthreshold',   ''); % see also below
cfg.scalpthreshold   = ft_getopt(cfg, 'scalpthreshold',   ''); % see also below
% earlier version of smooth and threshold specification
cfg.smooth      = ft_getopt(cfg, 'smooth',      '');
cfg.threshold   = ft_getopt(cfg, 'threshold',   '');

% chech whether earlier version of smooth and threshold was specified
if ~(isempty(cfg.smooth))
    if isempty(cfg.brainsmooth)
        cfg.brainsmooth=cfg.smooth;
        warning('Smoothing can be specified separately for scalp and brain. User-specified smoothing will be applied for brainmask.')
    end
    if isempty(cfg.scalpsmooth)
        cfg.scalpsmooth=cfg.smooth;
        warning('Smoothing can be specified separately for scalp and brain. User-specified smoothing will be applied for scalpmask.')
    end
end
if ~(isempty(cfg.threshold))
    if isempty(cfg.brainthreshold)
        cfg.brainthreshold=cfg.threshold;
        warning('Threshold can be specified separately for scalp and brain. User-specified threshold will be applied for brainmask.')
    end
    if isempty(cfg.scalpthreshold)
        cfg.scalpthreshold=cfg.threshold;
        warning('Threshold can be specified separately for scalp and brain. User-specified threshold  will be applied for scalpmask.')
        
    end
end
% then set defaults again
cfg.brainsmooth      = ft_getopt(cfg, 'brainsmooth',      5);
cfg.scalpsmooth      = ft_getopt(cfg, 'scalpsmooth',      5);
cfg.brainthreshold   = ft_getopt(cfg, 'brainthreshold',   0.5);
cfg.scalpthreshold   = ft_getopt(cfg, 'scalpthreshold',   0.1);

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

if ~isfield(cfg,'name') 
  if ~strcmp(cfg.write,'yes')
    tmp = tempname;
    cfg.name = tmp;
  else
    error('you must specify the output filename in cfg.name');
  end
end 

if ~iscell(cfg.output)
  % ensure it to be cell, to allow for multiple outputs
  cfg.output = {cfg.output};
end

% check whether spm is needed to generate tissue probability maps
if numel(cfg.output) == 1 && strcmp('scalp', cfg.output)
    needtpm= 0; % no tpm is necessary for creating (cummulative type) scalpmask
else
needtpm    = any(ismember(cfg.output, {'tpm' 'brain' 'skull' 'scalp' 'skullstrip'}));
end
hastpm     = isfield(mri, 'gray') && isfield(mri, 'white') && isfield(mri, 'csf');

if needtpm && ~hastpm
  % spm needs to be used for the creation of the tissue probability maps  
  dotpm = 1;
else
  dotpm = 0;
end

needana    = any(ismember(cfg.output, {'scalp' 'skullstrip'})) || dotpm;
hasanatomy = isfield(mri, 'anatomy');
if needana && ~hasanatomy
  error('the input volume needs an anatomy-field');
end

% perform optional downsampling before segmentation
if cfg.downsample ~= 1
  tmpcfg            = [];
  tmpcfg.downsample = cfg.downsample;
  tmpcfg.smooth     = 'no'; % smoothing is done in ft_volumesegment itself
  mri = ft_volumedownsample(tmpcfg, mri);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the tissue probability maps if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dotpm 
  
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
  
  Va = ft_write_mri([cfg.name,'.img'], mri.anatomy, 'transform', mri.transform, 'spmversion', cfg.spmversion);

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
        % put the transformation-matrix in the headers
        V(j).mat = mri.transform;
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
        % put the transformation-matrix in the headers
        V(j).mat = mri.transform;
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
  segment.anatomy   = mri.anatomy;

  % flip the volumes back according to the changes introduced by align_ijk2xyz
  for k = 1:3
    if flipflags(k)
      segment.gray    = flipdim(segment.gray, k);
      segment.anatomy = flipdim(segment.anatomy, k);
      if isfield(segment, 'white'), segment.white = flipdim(segment.white, k); end
      if isfield(segment, 'csf'),   segment.csf   = flipdim(segment.csf, k);   end
    end
  end

  if ~all(permutevec == [1 2 3])
    segment.gray    = ipermute(segment.gray,    permutevec);
    segment.anatomy = ipermute(segment.anatomy, permutevec);
    if isfield(segment, 'white'), segment.white = ipermute(segment.white, permutevec); end
    if isfield(segment, 'csf'),   segment.csf   = ipermute(segment.csf,   permutevec); end
    segment.dim  = size(segment.gray);
  end

else
  % rename the data
  segment = mri;
  clear mri;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now the data contains the tissue probability maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the requested output fields


removefields = {'anatomy' 'csf' 'gray' 'white'};

% check if smoothing or thresholding is required

dosmooth_brain = ~strcmp(cfg.brainsmooth,'no');
if dosmooth_brain && ischar(cfg.brainsmooth)
    error('invalid value %s for cfg.brainsmooth', cfg.brainsmooth);
end
dosmooth_scalp = ~strcmp(cfg.scalpsmooth, 'no');
if dosmooth_scalp && ischar(cfg.scalpsmooth)
    error('invalid value %s for cfg.scalpsmooth', cfg.scalpsmooth);
end

dothres_brain = ~strcmp(cfg.brainthreshold,'no');
if dothres_brain && ischar(cfg.brainthreshold)
    error('invalid value %s for cfg.brainthreshold', cfg.brainthreshold);
end
dothres_scalp = ~strcmp(cfg.scalpthreshold,'no');
if dothres_scalp && ischar(cfg.scalpthreshold)
    error('invalid value %s for cfg.scalpthreshold', cfg.scalpthreshold);
end
 
outp = cfg.output;

if any(strcmp(outp, 'tpm'))
    % output: probability tissue maps
    removefields = intersect(removefields, {'anatomy'});   
end

if any(strcmp(outp, 'brain')) || any(strcmp(outp, 'skull')) || any(strcmp(outp,'scalp')) || any(strcmp(outp,'skullstrip'))
    
    createoutputs=true;
    while createoutputs
        %create scalpmask - no tpm or brainmask is required to create it
        if any(strcmp('scalp',outp))
            
            fprintf('creating scalpmask\n');
            anatomy = segment.anatomy;
            if dosmooth_scalp, anatomy = volumesmooth(anatomy,     cfg.scalpsmooth,    'anatomy');
            else fprintf('no smoothing applied on anatomy for scalp segmentation\n');
            end
            if dothres_scalp, anatomy = volumethreshold(anatomy,  cfg.scalpthreshold, 'anatomy');
            else fprintf('no threshold applied on anatomy for scalp segmentation\n')
            end
            anatomy = volumefillholes(anatomy, 2); % FIXME why along the second dimension?
            
            scalpmask = anatomy>0;
            clear anatomy;
            
            % output: scalp (cummulative) (if this is the only requested
            % output)
            if numel(outp) ==1
                segment.scalp = scalpmask;
                
                break
            end
        end
            
        % create the brain from the tpm
        fprintf('creating brainmask\n');
        brain = segment.gray + segment.white + segment.csf;
        if dosmooth_brain, brain = volumesmooth(brain,  cfg.brainsmooth, 'brainmask');
        else fprintf('no smoothing applied on brainmask\n')
        end
        if dothres_brain, brain = volumethreshold(brain, cfg.brainthreshold, 'brainmask');
        else fprintf('no threshold applied on brainmask\n')
        end
        
        % output: skullstrip
        if any(strcmp('skullstrip',outp))
            
            
            fprintf('creating skullstripped anatomy\n');
            brain_ss = cast(brain, class(segment.anatomy));
            segment.anatomy = segment.anatomy.*brain_ss;
            clear brain_ss;
            removefields = intersect(removefields, {'gray' 'white' 'csf'});
            if numel(outp) == 1
               
                break
            end
        end
        % make binary mask from brain
        brainmask = brain>0;
        clear brain;
        
        
        % output: brain
        if any(strcmp(outp,'brain'))
            segment.brain = brainmask;
            
            if numel(outp) ==1
                break
            end
        end
        
        % output: skull
        if any(strcmp('skull', outp)) || any(strcmp('scalp', outp))
            % create skull from brain mask FIXME check this (e.g. strel_bol)
            fprintf('creating skullmask\n');
            braindil      = imdilate(brainmask>0, strel_bol(6));
            skullmask = braindil & ~brainmask;
            if any(strcmp(outp, 'skull'))
                segment.skull = skullmask;
                if numel(outp) == 1
                    break
                end
            end
            clear braindil;
            
            
            % output: scalp (exclusive type)
            if numel(outp) > 1 && any(strcmp('scalp',outp))
                scalpmask(brainmask>0)=0;
                clear brainmask;
                scalpmask(skullmask>0)=0;
                clear skullmask;
                segment.scalp=scalpmask;
                clear scalpmask;
            end
        end
        
        createoutputs=false; % exit while loop
    end
elseif strcmp(outp, 'tpm')
    %
else
    error('unknown output %s requested', cfg.output);
end

% remove unnecessary fields
for k = 1:numel(removefields)
  if isfield(segment, removefields{k})
    segment = rmfield(segment, removefields{k}); 
  end
end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble trackconfig
ft_postamble callinfo
ft_postamble previous mri
ft_postamble history segment
ft_postamble savevar segment

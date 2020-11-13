function [object] = ft_convert_coordsys(object, target, varargin)

% FT_CONVERT_COORDSYS changes the coordinate system of the input object to the
% specified coordinate system. The coordinate system of the input object is
% determined from the 'coordsys' field in the input data, or needs to be determined
% and specified interactively by the user.
%
% Use as
%   [output] = ft_convert_coordsys(input, target)
%   [output] = ft_convert_coordsys(input, target, method)
%   [output] = ft_convert_coordsys(input, target, method, template)
% to determine and convert the coordinate system.
%
% With the optional method input argument you can determine whether to use SPM for an
% affine or non-linear transformation.
%   method = 0: only an approximate coregistration (default for non-MRI data)
%   method = 1: an approximate coregistration, followed by spm_affreg
%   method = 2: an approximate coregistration, followed by spm_normalise (default for MRI data)
%
% The following input data structures are supported
%   electrode or gradiometer array, see FT_DATATYPE_SENS
%   volume conduction model, see FT_DATATYPE_HEADMODEL
%   source model, see FT_DATATYPE_SOURCE and FT_PREPARE_SOURCEMODEL
%   anatomical mri, see FT_DATATYPE_VOLUME
%   segmented mri, see FT_DATATYPE_SEGMENTATION
%   anatomical or functional atlas, see FT_READ_ATLAS
%
% Possible input coordinate systems are 'ctf', 'bti', '4d', 'neuromag' and 'itab'.
% Possible target coordinate systems are 'acpc'.
%
% Note that the conversion will be an automatic and approximate conversion, not
% taking into account differences in individual anatomies/differences in conventions
% where to put the fiducials.
%
% See also FT_DETERMINE_COORDSYS, FT_DETERMINE_UNITS, FT_CONVERT_UNITS, FT_PLOT_AXES, FT_PLOT_XXX

% Undocumented options
%   feedback  = string, 'yes' or 'no' (default = 'no')

% Copyright (C) 2005-2020, Robert Oostenveld & Jan-Mathijs Schoffelen
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

if nargin>3 && isnumeric(varargin{1})
  % old-style with 2 extra input arguments
  tmp = {'method', varargin{1}, 'template', varargin{2}};
  varargin = tmp;
elseif nargin>2 && isnumeric(varargin{1})
  % old-style with 1 extra input argument
  tmp = {'method', varargin{1}};
  varargin = tmp;
end

method        = ft_getopt(varargin, 'method');    % default is handled below
templatefile  = ft_getopt(varargin, 'template');  % default is handled in the SPM section
feedback      = ft_getopt(varargin, 'feedback', 'no');

if isempty(method)
  if isfield(object, 'transform') && isfield(object, 'anatomy')
    % the default for an anatomical MRI is to start with an approximate alignment,
    % followed by a call to spm_normalise for a better quality alignment
    method = 2;
  else
    % the default for all other objects is to do only an approximate alignment
    method = 0;
  end
end

if isdeployed && method>0 && isempty(templatefile)
  ft_error('you need to specify a template filename for the coregistration');
end

if ~isfield(object, 'coordsys') || isempty(object.coordsys)
  % determine the coordinate system of the input object
  object = ft_determine_coordsys(object, 'interactive', 'yes');
  if ~isfield(object, 'coordsys') || isempty(object.coordsys)
    % the call to ft_determine_coordsys should have taken care of this, but
    % it is possible that the user aborted the coordinate system
    % determination. See http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2526
    ft_error('the coordinate system of the geometrical object is not specified');
  end
  ft_notice('the coordinate system is determined as ''%s''', object.coordsys);
end

if ~isfield(object, 'unit') || isempty(object.unit)
  % determine the units of the input object
  object = ft_determine_units(object);
  ft_notice('the units are determined as ''%s''', object.unit);
end

% all of the internal logic inside this function requires that the units are in millimeter
original = object;
object = ft_convert_units(object, 'mm');

if ~ismember(object.coordsys, {'spm', 'mni', 'fsaverage', 'tal'}) && ismember(target, {'spm', 'mni', 'fsaverage', 'tal'})
  % the input appears to be an individual subject MRI which has not been rescaled
  % the target is a template coordinate system
  % see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3304 and http://bit.ly/2sw7eC4
  ft_warning('Not applying any scaling, using ''acpc'' instead of ''%s''. See http://bit.ly/2sw7eC4', target);
  target = 'acpc';
end

% these are the 48 generic axis orientation triplets
%   a = anterior
%   p = posterior
%   l = left
%   r = right
%   s = superior
%   i = inferior

generic = {
  'als'; 'ali'; 'ars'; 'ari';...
  'pls'; 'pli'; 'prs'; 'pri';...
  'las'; 'lai'; 'ras'; 'rai';...
  'lps'; 'lpi'; 'rps'; 'rpi';...
  'asl'; 'ail'; 'asr'; 'air';...
  'psl'; 'pil'; 'psr'; 'pir';...
  'sal'; 'ial'; 'sar'; 'iar';...
  'spl'; 'ipl'; 'spr'; 'ipr';...
  'sla'; 'ila'; 'sra'; 'ira';...
  'slp'; 'ilp'; 'srp'; 'irp';...
  'lsa'; 'lia'; 'rsa'; 'ria';...
  'lsp'; 'lip'; 'rsp'; 'rip'}';

specific = {'ctf', 'bti', 'fourd', 'neuromag', 'itab', 'acpc', 'mni', 'spm', 'fsaverage', 'tal'};

% generic orientation triplets (like RAS and ALS) are not specific with regard to the origin
if ismember(object.coordsys, generic) && strcmp(target, 'acpc')
  if method==0
    ft_error('approximately converting from %s to %s is not supported', object.coordsys, target);
  elseif method>0
    % converting an anatomical MRI from RAS, ALS etc. to ACPC using SPM might work
  end
elseif ismember(object.coordsys, generic) && ~ismember(target, generic)
  % other conversions from generic orientation triplets (like RAS and ALS) are also not supported
  ft_error('converting from %s to %s is not supported', object.coordsys, target);
end


%--------------------------------------------------------------------------
% Start with an approximate alignment, this is based on transformation matrices
% that were determined by clicking on the CTF fiducial locations in the canonical
% T1 template MRI.
%
% All of the transformation matrices here are expressed in millimeter.

% this is based on the ear canals, see ALIGN_CTF2ACPC
acpc2ctf = [
  0.0000  0.9987  0.0517  34.7467
 -1.0000  0.0000  0.0000   0.0000
  0.0000 -0.0517  0.9987  52.2749
  0.0000  0.0000  0.0000   1.0000
  ];

% this is based on the ear canals, see ALIGN_NEUROMAG2ACPC
acpc2neuromag = [
  1.0000  0.0000  0.0000   0.0000
  0.0000  0.9987  0.0517  34.7467
  0.0000 -0.0517  0.9987  52.2749
  0.0000  0.0000  0.0000   1.0000
  ];

% see http://freesurfer.net/fswiki/CoordinateSystems
fsaverage2mni = [
  0.9975   -0.0073    0.0176   -0.0429
  0.0146    1.0009   -0.0024    1.5496
 -0.0130   -0.0093    0.9971    1.1840
  0.0000    0.0000    0.0000    1.0000
  ];

% this is a 90 degree rotation around the z-axis
ctf2neuromag = [
  0.0000   -1.0000    0.0000    0.0000
  1.0000    0.0000    0.0000    0.0000
  0.0000    0.0000    1.0000    0.0000
  0.0000    0.0000    0.0000    1.0000
  ];

% these are all combinations of 90 degree rotations and/or flips along one of the axes
for i=1:length(generic)
  for j=1:length(generic)
    xxx = generic{i};
    yyy = generic{j};
    eval(sprintf('%s2%s = transform_generic(''%s'', ''%s'');', yyy, xxx, xxx, yyy));
  end
end

% affine transformation from MNI to Talairach, see http://imaging.mrc-cbu.cam.ac.uk/imaging/MniTalairach
% the non-linear (i.e. piecewise linear) transform between MNI and Talairach are implemented elsewhere, see the functions MNI2TAL and TAL2MNI
mni2tal = [
  0.8800    0.0000    0.0000   -0.8000
  0.0000    0.9700    0.0000   -3.3200
  0.0000    0.0500    0.8800   -0.4400
  0.0000    0.0000    0.0000    1.0000
  ];

% the CTF and BTI coordinate system are the same, see http://www.fieldtriptoolbox.org/faq/how_are_the_different_head_and_mri_coordinate_systems_defined/
ctf2bti = eye(4);

% the Neuromag and Itab coordinate system are the same, see http://www.fieldtriptoolbox.org/faq/how_are_the_different_head_and_mri_coordinate_systems_defined/#details-of-the-ctf-coordinate-system
neuromag2itab = eye(4);

% BTI and 4D are different names for the same system
bti2fourd = eye(4);

% the SPM and MNI coordinate system are the same, see http://www.fieldtriptoolbox.org/faq/acpc/
spm2mni = eye(4);

% the SPM (aka MNI) and ACPC coordinate system are not the same but similar enough, see http://www.fieldtriptoolbox.org/faq/acpc/
spm2acpc = eye(4);
mni2acpc = eye(4);

% the CTF, BTI and 4D coordinate systems are all ALS coordinate systems
% but the origin is poorly defined in ALS, hence converting from ALS to another is problematic
ctf2als   = eye(4);
bti2als   = eye(4);
fourd2als = eye(4);

% the Neuromag, Itab, ACPC, MNI, SPM and FSAVERAGE coordinate systems are all RAS coordinate systems
% but the origin is poorly defined in RAS, hence converting from RAS to another is problematic
neuromag2ras  = eye(4);
itab2ras      = eye(4);
acpc2ras      = eye(4);
mni2ras       = eye(4);
spm2ras       = eye(4);
fsaverage2ras = eye(4);
tal2ras       = eye(4);

% make the combined and the inverse transformations where possible
coordsys = [specific generic];
implemented = zeros(length(coordsys)); % this is only for debugging
for i=1:numel(coordsys)
  for j=1:numel(coordsys)
    xxx = coordsys{i};
    yyy = coordsys{j};
    
    if isequal(xxx, yyy)
      % construct the transformations on the diagonal
      eval(sprintf('%s2%s = eye(4);', xxx, yyy));
      implemented(i,j) = 1;
    elseif exist(sprintf('%s2%s', xxx, yyy), 'var')
      % construct the inverse transformations
      eval(sprintf('%s2%s = inv(%s2%s);', yyy, xxx, xxx, yyy));
      implemented(i,j) = 2;
      implemented(j,i) = 2;
    elseif ismember(xxx, specific) && ismember(yyy, generic)
      % try to make the transformation (and inverse) with a two-step approach
      % since we go from specific to generic and thereby loose the origin information anyway, it is fine to use any intermediate step
      for k=1:numel(coordsys)
        zzz = coordsys{k};
        if exist(sprintf('%s2%s', xxx, zzz), 'var') && exist(sprintf('%s2%s', zzz, yyy), 'var')
          eval(sprintf('%s2%s = %s2%s * %s2%s;', xxx, yyy, zzz, yyy, xxx, zzz));
          eval(sprintf('%s2%s = inv(%s2%s);', yyy, xxx, xxx, yyy));
          implemented(i,j) = 3;
          implemented(j,i) = 3;
          break
        end
      end % for k
    elseif ismember(xxx, specific) && ismember(yyy, specific)
      % try to make the transformation (and inverse) with a two-step approach
      % do not use the generic orientation triplets (like RAS and ALS) as intermediate steps between two specific coordinate systems
      for k=1:numel(specific)
        zzz = specific{k};
        if exist(sprintf('%s2%s', xxx, zzz), 'var') && exist(sprintf('%s2%s', zzz, yyy), 'var')
          eval(sprintf('%s2%s = %s2%s * %s2%s;', xxx, yyy, zzz, yyy, xxx, zzz));
          eval(sprintf('%s2%s = inv(%s2%s);', yyy, xxx, xxx, yyy));
          implemented(i,j) = 3;
          implemented(j,i) = 3;
          break
        end
      end % for k
    end
    
  end % for j
end % for i

% these conversions should be done using FT_VOLUMENORMALISE, as they imply scaling
clear acpc2spm acpc2mni acpc2fsaverage acpc2tal

% converting to/from TAL is only possible for some specific template coordinate systems
clear bti2tal ctf2tal fourd2tal itab2tal neuromag2tal
clear tal2bti tal2ctf tal2fourd tal2itab tal2neuromag

% the origin is poorly defined in generic orientation triplets (like RAS and ALS), hence converting them to any specific coordinate system is problematic
% the only conversion supported is from generic orientation triplets (like RAS and ALS) to ACPC, and only when using SPM
if method<1
  for i=1:length(generic)
    for j=1:length(specific)
      xxx = generic{i};
      yyy = specific{j};
      eval(sprintf('clear %s2%s', xxx, yyy));
    end
  end
end

% this is only for debugging the coverage of conversions, note that some of them are deleted further down
if false
  % update the list of implemented transformations, since some might have been cleared
  for i=1:length(coordsys)
    for j=1:length(coordsys)
      xxx = coordsys{i};
      yyy = coordsys{j};
      if ~exist(sprintf('%s2%s', xxx, yyy), 'var')
        implemented(i,j) = 0;
      end
    end
  end
  figure; imagesc(implemented); caxis([0 3]);
  xticklabels(coordsys); xticks(1:numel(coordsys));
  yticklabels(coordsys); yticks(1:numel(coordsys));
end

if strcmp(object.coordsys, '4d')
  xxx = 'fourd'; % '4d' is not a valid variable name
else
  xxx = object.coordsys;
end

if strcmp(target, '4d')
  yyy = 'fourd'; % '4d' is not a valid variable name
else
  yyy = target;
end

if exist(sprintf('%s2%s', xxx, yyy), 'var')
  transform = eval(sprintf('%s2%s', xxx, yyy));
  object = ft_transform_geometry(transform, object);
  object.coordsys = target;
else
  ft_error('converting from %s to %s is not supported', object.coordsys, target);
end


%--------------------------------------------------------------------------
% Do a second round of affine registration (rigid body) to get improved
% alignment with ACPC coordinate system. This is needed because there may be
% different conventions defining LPA and RPA. The affine registration may
% fail however, e.g. if the initial alignment is not close enough. In that
% case SPM will throw an error.
%
% We expect the template MRIs to be expressed in millimeter and to be
% approximately aligned with ACPC.

if method>0
  if ~isfield(object, 'transform') || ~isfield(object, 'anatomy')
    ft_error('affine or non-linear transformation are only supported for anatomical MRIs');
  end
  if ~isfield(object, 'unit') || ~strcmp(object.unit, 'mm')
    ft_error('affine or non-linear transformation require the anatomial MRI to be expressed in mm');
  end
  if ~any(ismember(object.coordsys, {'acpc', 'spm', 'mni', 'fsaverage', 'tal'}))
    % this constraint could be relaxed if we would know that the template is expressed in another coordinate system
%     ft_error('affine or non-linear transformation is only supported for data in an SPM-like coordinate systems');
  end
  
  % this requires SPM to be on the path. However, this is not the proper place to
  % choose between SPM versions. The user can either use cfg.spmversion in a high-level
  % function, or has to add the path to the desired SPM version by hand.
  ft_hastoolbox('spm', -1);
end


if method==1
  % use spm_affreg
  if isempty(templatefile)
    switch lower(spm('ver'))
      case 'spm2'
        if isdeployed
          if nargin<3, ft_error('you need to specify a template filename when in deployed mode and using method==1'); end
        else
          templatefile = fullfile(spm('Dir'),'templates','T1.mnc');
        end
      case 'spm8'
        if isdeployed
          if nargin<3, ft_error('you need to specify a template filename when in deployed mode and using method==1'); end
        elseif isempty(templatefile)
          templatefile = fullfile(spm('Dir'),'templates','T1.nii');
        end
      case 'spm12'
        if isdeployed
          if nargin<3, ft_error('you need to specify a template filename when in deployed mode and using method==1'); end
        else
          templatefile = fullfile(spm('Dir'),'toolbox','OldNorm','T1.nii');
          if ~exist('spm_affreg', 'file')
            addpath(fullfile(spm('Dir'),'toolbox','OldNorm'));
          end
        end
        fprintf('using ''OldNorm'' affine registration\n');
      otherwise
        ft_error('unsupported SPM version');
    end
  end
  template = ft_read_mri(templatefile);
  
  tname1 = [tempname, '.img'];
  tname2 = [tempname, '.img'];
  V1 = ft_write_mri(tname1, object.anatomy,  'transform', object.transform,  'spmversion', spm('ver'), 'dataformat', 'nifti_spm');
  V2 = ft_write_mri(tname2, template.anatomy, 'transform', template.transform, 'spmversion', spm('ver'), 'dataformat', 'nifti_spm');
  
  % the below, using just spm_affreg does not work robustly enough in some cases
  flags.regtype = 'rigid';
  [M, scale]    = spm_affreg(V1,V2,flags);
  
  % some juggling around with the transformation matrices
  vox2head    = object.transform;
  vox2acpc2  = M \ V1.mat;
  acpc2head2 = vox2head / vox2acpc2;
  
  % update the transformation matrix
  object.transform = vox2acpc2;
  
  % delete the temporary files
  delete(tname1); delete(strrep(tname1, 'img', 'hdr'));
  delete(tname2); delete(strrep(tname2, 'img', 'hdr'));
  
elseif method==2
  % use spm_normalise
  if isempty(templatefile)
    switch lower(spm('ver'))
      case 'spm2'
        if isdeployed
          if nargin<3, ft_error('you need to specify a template filename when in deployed mode and using method==2'); end
        else
          templatefile = fullfile(spm('Dir'),'templates','T1.mnc');
        end
      case 'spm8'
        if isdeployed
          if nargin<3, ft_error('you need to specify a template filename when in deployed mode and using method==2'); end
        else
          templatefile = fullfile(spm('Dir'),'templates','T1.nii');
        end
      case 'spm12'
        % this uses the 'OldNorm' functionality, so the path needs to be added, can only be done if non-deployed.
        if isdeployed
          if nargin<3, ft_error('you need to specify a template filename when in deployed mode and using method==2'); end
        else
          templatefile = fullfile(spm('Dir'),'toolbox','OldNorm','T1.nii');
          if ~exist('spm_normalise', 'file')
            addpath(fullfile(spm('Dir'),'toolbox','OldNorm'));
          end
        end
        fprintf('using ''OldNorm'' normalisation\n');
      otherwise
          ft_error('unsupported SPM version');
    end
  end
  template = ft_read_mri(templatefile);
  
  tname1 = [tempname, '.img'];
  tname2 = [tempname, '.img'];
  V1 = ft_write_mri(tname1, object.anatomy,  'transform', object.transform,  'spmversion', spm('ver'), 'dataformat', 'nifti_spm');
  V2 = ft_write_mri(tname2, template.anatomy, 'transform', template.transform, 'spmversion', spm('ver'), 'dataformat', 'nifti_spm');
  
  flags.nits       = 0; % set number of non-linear iterations to zero
  flags.regtype    = 'rigid';
  params           = spm_normalise(V2,V1,[],[],[],flags);
  
  % some juggling around with the transformation matrices
  vox2head   = object.transform;
  acpc2head2 = V1.mat*params.Affine/V2.mat;
  vox2acpc2  = acpc2head2\vox2head;
  
  % update the transformation matrix
  object.transform = vox2acpc2;
  
  % delete the temporary files
  delete(tname1); delete(strrep(tname1, 'img', 'hdr'));
  delete(tname2); delete(strrep(tname2, 'img', 'hdr'));
end

if istrue(feedback)
  % give some graphical feedback
  ft_determine_coordsys(object, 'interactive', 'no', 'fontsize', 15);
  % also add the original axes
  if method==0
    ft_plot_axes([], 'transform', transform, 'unit', 'mm', 'coordsys', original.coordsys, 'fontsize', 15);
  elseif method>0
    transform = object.transform * inv(original.transform);
    ft_plot_axes([], 'transform', transform, 'unit', 'mm', 'coordsys', original.coordsys, 'fontsize', 15);
  end
end

% all of the internal logic inside this function requires that the units are in millimeter
% convert back to the original units
object = ft_convert_units(object, original.unit);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION to construct generic transformations such as RAS2ALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = transform_generic(from, to)

ap_in  = find(from=='a' | from=='p');
ap_out = find(to=='a'   | to=='p');
lr_in  = find(from=='l' | from=='r');
lr_out = find(to=='l'   | to=='r');
si_in  = find(from=='s' | from=='i');
si_out = find(to=='s'   | to=='i');

% index axis according to ap,lr,si
order_in  = [ap_in  lr_in  si_in];
order_out = [ap_out lr_out si_out];

% check whether one of the axis needs flipping
flip = 2.*(0.5-double(from(order_in)~=to(order_out)));

T = zeros(4);
for k = 1:3
  T(order_out(k),order_in(k)) = flip(k);
end
T(4,4) = 1;

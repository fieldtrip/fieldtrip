function [object] = ft_convert_coordsys(object, target, method, templatefile)

% FT_CONVERT_COORDSYS changes the coordinate system of the input object to
% the specified coordinate system. The coordinate system of the input
% object is determined from the structure field object.coordsys, or needs to
% be determined and specified interactively by the user.
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

if ~isfield(object, 'coordsys') || isempty(object.coordsys)
  % determine the coordinate system of the input object
  object = ft_determine_coordsys(object, 'interactive', 'yes');
end

if ~isfield(object, 'coordsys') || isempty(object.coordsys)
  % the call to ft_determine_coordsys should have taken care of this, but
  % it is possible that the user aborted the coordinate system
  % determination. See http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2526
  ft_error('the coordinate system of the geometrical object is not specified');
end

if any(strcmp(target, {'spm', 'mni', 'tal'}))
  % see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3304
  ft_warning('Not applying any scaling, using ''acpc'' instead of ''%s''. See http://bit.ly/2sw7eC4', target);
  target = 'acpc';
end

if nargin<3
  if isfield(object, 'transform') && ~isfield(object, 'anatomy')
    % the default for an MRI is to start with an approximate alignment,
    % followed by a call to spm_normalise for a better quality alignment
    method = 2;
  else
    % the default for all other objects is to do only an approximate alignment
    method = 0;
  end
end

if isdeployed && method>0
  needtemplate = true;
else
  needtemplate = false;
end

hastemplate = nargin>3;
if needtemplate && ~hastemplate
  ft_error('you need to specify a template filename for the coregistration');
end

%--------------------------------------------------------------------------
% Start with an approximate alignment, this is based on transformation matrices
% that were determined by clicking on the CTF fiducial locations in the canonical
% T1 template MRI

if ~strcmpi(target, object.coordsys)
  
  acpc2ctf = [
    0.0000  0.9987  0.0517  34.7467
    -1.0000  0.0000  0.0000   0.0000
    0.0000 -0.0517  0.9987  52.2749
    0.0000  0.0000  0.0000   1.0000
    ];
  
  acpc2neuromag = [
    1.0000  0.0000  0.0000   0.0000
    0.0000  0.9987  0.0517  34.7467
    0.0000 -0.0517  0.9987  52.2749
    0.0000  0.0000  0.0000   1.0000
    ];
  
  % these are alternative names
  acpc2itab   = acpc2neuromag;
  acpc2bti    = acpc2ctf;
  acpc2fourd  = acpc2ctf;
  
  % also allow reverse coordinate system conversions
  ctf2acpc      = inv(acpc2ctf);
  neuromag2acpc = inv(acpc2neuromag);
  itab2acpc     = inv(acpc2itab);
  bti2acpc      = inv(acpc2bti);
  fourd2acpc    = inv(acpc2fourd);
  
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
    ft_error('conversion from %s to %s is not supported', object.coordsys, target);
  end
  
end % approximate alignment

%--------------------------------------------------------------------------
% Do a second round of affine registration (rigid body) to get improved
% alignment with ACPC coordinate system. this is needed because there may be
% different conventions defining LPA and RPA. The affine registration may
% fail however, e.g. if the initial alignment is not close enough. In that
% case SPM will throw an error

if method>0
  if ~isfield(object, 'transform') || ~isfield(object, 'anatomy')
    ft_error('affine or non-linear transformation using SPM are only supported for anatomical MRIs');
  end
  if ~isfield(object, 'unit') || ~strcmp(object.unit, 'mm')
    ft_error('affine or non-linear transformation using SPM require the anatomial MRI to be expressed in mm');
  end
  if ~strcmp(object.coordsys, 'acpc')
    % this constraint could be relaxed if we would know that the template is expressed in another coordinate system
    ft_error('affine or non-linear transformation using SPM are only supported for ACPC');
  end
  
  % this requires SPM to be on the path. However, this is not the proper place to
  % choose between SPM versions. The user can either use cfg.spmversion in a high-level
  % function, or has to add the path to the desired SPM version by hand.
  ft_hastoolbox('spm', -1);
end


if method==1
  % use spm_affreg
  
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
      else
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
  template = ft_read_mri(templatefile);
  
  tname1 = [tempname, '.img'];
  tname2 = [tempname, '.img'];
  V1 = ft_write_mri(tname1, object.anatomy,  'transform', object.transform,  'spmversion', spm('ver'), 'dataformat', 'nifti_spm');
  V2 = ft_write_mri(tname2, template.anatomy, 'transform', template.transform, 'spmversion', spm('ver'), 'dataformat', 'nifti_spm');
  
  % the below, using just spm_affreg does not work robustly enough in some cases
  flags.regtype = 'rigid';
  [M, scale]    = spm_affreg(V1,V2,flags);
  
  % some juggling around with the transformation matrices
  mrivox2mrihead    = object.transform;
  mrivox2acpchead2  = M \ V1.mat;
  acpchead2mrihead2 = mrivox2mrihead / mrivox2acpchead2;
  
  % update the transformation matrix
  object.transform     = mrivox2acpchead2;
  
  % this one is unchanged
  object.vox2headOrig  = mrivox2mrihead;
  
  % these are new
  object.vox2head      = mrivox2acpchead2;
  object.head2headOrig = acpchead2mrihead2;
  
  % delete the temporary files
  delete(tname1); delete(strrep(tname1, 'img', 'hdr'));
  delete(tname2); delete(strrep(tname2, 'img', 'hdr'));
  
elseif method==2
  % use spm_normalise
  
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
  template = ft_read_mri(templatefile);
  
  tname1 = [tempname, '.img'];
  tname2 = [tempname, '.img'];
  V1 = ft_write_mri(tname1, object.anatomy,  'transform', object.transform,  'spmversion', spm('ver'), 'dataformat', 'nifti_spm');
  V2 = ft_write_mri(tname2, template.anatomy, 'transform', template.transform, 'spmversion', spm('ver'), 'dataformat', 'nifti_spm');
  
  flags.nits       = 0; %set number of non-linear iterations to zero
  flags.regtype    = 'rigid';
  params           = spm_normalise(V2,V1,[],[],[],flags);
  
  % some juggling around with the transformation matrices
  mrivox2mrihead    = object.transform;
  acpchead2mrihead2 = V1.mat*params.Affine/V2.mat;
  mrivox2acpchead2  = acpchead2mrihead2\mrivox2mrihead;
  
  % update the transformation matrix
  object.transform     = mrivox2acpchead2;
  
  % this one is unchanged
  object.vox2headOrig  = mrivox2mrihead;
  
  % these are new
  object.vox2head      = mrivox2acpchead2;
  object.head2headOrig = acpchead2mrihead2;
  
  % delete the temporary files
  delete(tname1); delete(strrep(tname1, 'img', 'hdr'));
  delete(tname2); delete(strrep(tname2, 'img', 'hdr'));
end

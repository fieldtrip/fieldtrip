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
% Recognized and supported coordinate systems are 'ctf', 'bti', '4d', 'yokogawa',
% 'eeglab', 'neuromag', 'itab', 'acpc', 'spm', 'mni', 'fsaverage', 'tal', 'scanras',
% 'scanlps', 'dicom'.
%
% Furthermore, supported coordinate systems that do not specify the origin are 'ras',
% 'als', 'lps', etc. See https://www.fieldtriptoolbox.org/faq/coordsys for more
% details.
%
% Note that the conversion will be an automatic and approximate conversion, not
% taking into account differences in individual anatomies/differences in conventions
% where to put the fiducials.
%
% See also FT_DETERMINE_COORDSYS, FT_DETERMINE_UNITS, FT_CONVERT_UNITS, FT_PLOT_AXES, FT_PLOT_XXX

% Undocumented options
%   feedback  = string, 'yes' or 'no' (default = 'no')

% Copyright (C) 2005-2021, Robert Oostenveld & Jan-Mathijs Schoffelen
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

%--------------------------------------------------------------------------
% Do an initial affine registration (rigid body) alignment to the target coordinate
% system. This deals with 90-degree rotations between RAS and ALS, or between CTF and
% NEUROMAG, but also with approximate alignment between external NAS/LPA/RPA
% landmark-based coordinates to internal ACPCP landmark-based coordinate systems.

% these are the 48 generic axis orientation triplets, these specify the axes but no origin
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

if ismember(object.coordsys, generic) && strcmp(target, 'acpc') && method>0
  % converting anatomical MRI data from a generic coordinate system to ACPC with SPM should work
  % do an initial alignment of the anatomical MRI data to RAS to bring it closer to ACPC
  initial = ft_affinecoordinates(object.coordsys, 'ras');
else
  initial = ft_affinecoordinates(object.coordsys, target);
end

object = ft_transform_geometry(initial, object);
object.coordsys = target;

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

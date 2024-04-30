function mri = ft_defacevolume(cfg, mri)

% FT_DEFACEVOLUME allows you to de-identify an anatomical MRI by erasing specific
% regions, such as the face and ears. The interactive graphical user interface allows
% you to position a box over the anatomical data inside which all anatomical voxel
% values will be replaced by zero. You might have to call this function multiple
% times when both face and ears need to be removed. Following defacing, you should
% check the result with FT_SOURCEPLOT.
%
% Use as
%   mri = ft_defacevolume(cfg, mri)
%
% The configuration can contain the following options
%   cfg.method     = 'box', 'plane', 'spm' (default = 'box')
%
% If you specify the box method, the following options apply
%   cfg.translate  = initial position of the center of the box, or a point
%                    on the plane, (default = [0 0 0])
%   cfg.scale      = initial size of the box along each dimension (default is automatic)
%   cfg.rotate     = initial rotation of the box, or the plane (default = [0 0 0])
%   cfg.selection  = which voxels to keep, can be 'inside' or 'outside' (default = 'outside')
%   cfg.smooth     = 'no' or the FWHM of the gaussian kernel in voxels (default = 'no')
%   cfg.keepbrain  = 'no' or 'yes', segment and retain the brain (default = 'no')
%   cfg.feedback   = 'no' or 'yes', whether to provide graphical feedback (default = 'no')
%
% If you specify no smoothing, the selected area will be zero-masked. If you
% specify a certain amount of smoothing (in voxels FWHM), the selected area will
% be replaced by a smoothed version of the data.
%
% The spm method does not have any options, it uses SPM_DEFACE from the
% SPM12 toolbox.
%
% See also FT_ANONYMIZEDATA, FT_DEFACEMESH, FT_ANALYSISPIPELINE, FT_SOURCEPLOT

% Copyright (C) 2015-2022, Robert Oostenveld
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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar    mri
ft_preamble provenance mri

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% for backward compatibility
cfg = ft_checkconfig(cfg, 'renamedval', {'method', 'interactive', 'box'});

% set the defaults
cfg.method         = ft_getopt(cfg, 'method', 'box');
cfg.rotate         = ft_getopt(cfg, 'rotate', [0 0 0]);
cfg.scale          = ft_getopt(cfg, 'scale'); % the automatic default is determined further down
cfg.translate      = ft_getopt(cfg, 'translate', [0 0 0]);
cfg.transformorder = ft_getopt(cfg, 'transformorder', {'scale', 'rotate', 'translate'}); % T*R*S
cfg.selection      = ft_getopt(cfg, 'selection', 'outside');
cfg.smooth         = ft_getopt(cfg, 'smooth', 'no');
cfg.keepbrain      = ft_getopt(cfg, 'keepbrain', 'no');
cfg.feedback       = ft_getopt(cfg, 'feedback', 'no');

ismri  = ft_datatype(mri, 'volume') && isfield(mri, 'anatomy');
ismesh = isfield(mri, 'pos'); % triangles are optional

if ismri
  % check if the input data is valid for this function
  mri = ft_checkdata(mri, 'datatype', 'volume', 'feedback', 'yes');
end

switch cfg.method
  case 'spm'
    % this requires SPM12 on the path
    ft_hastoolbox('spm12', 1);

    % defacing relies on coregistration, which relies on the MRI being reasonably aligned for SPM
    mri = ft_checkdata(mri, 'hascoordsys', 'yes');

    % remember the original transformation matrix and coordinate system
    original = [];
    original.transform = mri.transform;
    original.coordsys  = mri.coordsys;
    mri = ft_convert_coordsys(mri, 'acpc');

    filename1 = {[tempname '.nii']};
    ft_write_mri(filename1{1}, mri, 'dataformat', 'nifti');

    % % apply a least squares pre-alignment step in order to make spm_deface more robust
    % % this could be done conditional on the modality/contrast, which is part of the BIDS filename
    % template = spm_vol(fullfile(spm('Dir'),'canonical','avg152PD.nii'));
    % template = spm_vol(fullfile(spm('Dir'),'canonical','avg152T1.nii'));
    % template = spm_vol(fullfile(spm('Dir'),'canonical','avg152T2.nii'));
    % filevol = spm_vol(filename1{1});
    % M = spm_affreg(template, filevol);
    % spm_get_space(filename1{1}, M * filevol.mat);

    filename2 = spm_deface(filename1);
    mri = ft_read_mri(filename2{1});

    % put the original transformation matrix and coordinate system back
    mri.transform = original.transform;
    mri.coordsys = original.coordsys;

    % clean up the temporary files
    delete(filename1{1});
    delete(filename2{1});

  case {'box' 'plane'}
    % this is an alternative implementation of the interactive method usinf FT_INTERACTIVEREALIGN
    % it aligns a box or a plane to the MRI or mesh, and then removes the points inside that box,
    % or below the plane

    if ismri
      % enhance the contrast of the volumetric data, see also FT_VOLUMEREALIGN
      dat  = double(mri.anatomy);
      dum  = unique(dat(:));
      dmin = dum(round(0.05*numel(dum))); % take the 5% value of the histogram
      dmax = dum(round(0.95*numel(dum))); % take the 95% value of the histogram
      dat  = (dat-dmin)./(dmax-dmin);
      mri.anatomy = dat;
    end

    if isequal(cfg.method, 'box')
      % construct a box with a unit length expressed in the units of the input mri or mesh
      % rather than using a triangulation, this specifies polygons for each of the 6 edges of the box
      box.unit = mri.unit;
      box.pos = [
        1  1  1
        1 -1  1
        -1 -1  1
        -1  1  1
        1  1 -1
        1 -1 -1
        -1 -1 -1
        -1  1 -1
        ]/2;
      box.poly = [
        1 2 3 4
        1 5 6 2
        2 6 7 3
        3 7 8 4
        4 8 5 1
        5 8 7 6
        ];

      % the default is to scale the box to 75 mm (or equivalent)
      defaultscale = [75 75 75] * ft_scalingfactor('mm', mri.unit);
    
    elseif isequal(cfg.method, 'plane')
      % call it box, make a plane
      box.unit = mri.unit;
      box.pos  = [ 1  1 0
                  -1  1 0
                  -1 -1 0
                   1 -1 0
                   0  0 0
                   0  0 -0.5];
      box.poly = [1 2 3 4];
      box.line = [5 6];

      % the default is to draw the plane as a 200x200 mm plane (or
      % equivalent, note that a deviation of this default does not have
      % functional consequences
      defaultscale = [200 200 1000] * ft_scalingfactor('mm', mri.unit);
    else
      ft_error('you can only specify a box or a plane as exclusion criterion');
    end

    surfaceonly = isfield(mri, 'tet') | isfield(mri, 'hex'); % only for tetrahedral or hexahedral meshes

    tmpcfg = keepfields(cfg, {'scale', 'rotate', 'translate', 'transformorder'});
    tmpcfg.scale = ft_getopt(cfg, 'scale', defaultscale);
    tmpcfg.showapply = 'no'; % do not show the apply button
    tmpcfg.template.axes = 'yes';
    if ismri
      tmpcfg.showlight = 'no';
      tmpcfg.showalpha = 'no';
      tmpcfg.template.mri = mri;
    elseif ismesh
      tmpcfg.showlight = 'yes';
      tmpcfg.showalpha = 'yes';
      tmpcfg.template.mesh = mri; % the input variable is called "mri" but it contains a mesh
      tmpcfg.template.meshstyle.facefolor = 'skin_medium';
      tmpcfg.template.meshstyle.edgecolor = 'none';
    end
    tmpcfg.individual.mesh = box;
    tmpcfg.individual.meshstyle = {'edgecolor', 'k', 'facecolor', 'y', 'facealpha', 0.3, 'surfaceonly', surfaceonly};
    tmpcfg = ft_interactiverealign(tmpcfg);

    % remember these for potential reuse outside of this function
    cfg.rotate    = tmpcfg.rotate;
    cfg.scale     = tmpcfg.scale;
    cfg.translate = tmpcfg.translate.*ft_scalingfactor('mm',mri.unit);

    % the template remains fixed, the individual is moved around
    R = rotate(cfg.rotate);
    T = translate(cfg.translate);
    if isequal(cfg.method, 'box')
      S = scale(cfg.scale);
    else isequal(cfg.method, 'plane')
      S = eye(4); % no scaling needs to be performed
    end
    % this is the transformation to get from the individual to the template
    transform = combine_transform(R, S, T, cfg.transformorder);

    if ismri
      % rather than converting the box to the MRI, do it the other way around
      [X, Y, Z] = ndgrid(1:mri.dim(1), 1:mri.dim(2), 1:mri.dim(3));
      voxpos = ft_warp_apply(mri.transform, [X(:) Y(:) Z(:)]);  % voxel positions in head coordinates
      voxpos = ft_warp_apply(inv(transform), voxpos);           % voxel positions in box coordinates

      remove = ...
        voxpos(:,1) > -0.5 & ...
        voxpos(:,1) < +0.5 & ...
        voxpos(:,2) > -0.5 & ...
        voxpos(:,2) < +0.5 & ...
        voxpos(:,3) > -0.5 & ...
        voxpos(:,3) < +0.5;

    elseif ismesh || issource
      % rather than converting the box to the mesh, do it the other way around
      meshpos = ft_warp_apply(inv(transform), mri.pos);         % vertex positions in box coordinates

      if isequal(cfg.method, 'box')
        remove = ...
          meshpos(:,1) > -0.5 & ...
          meshpos(:,1) < +0.5 & ...
          meshpos(:,2) > -0.5 & ...
          meshpos(:,2) < +0.5 & ...
          meshpos(:,3) > -0.5 & ...
          meshpos(:,3) < +0.5;
      else isequal(cfg.method, 'plane')
        remove = meshpos(:,3) < 0;
      end

    end

    if strcmp(cfg.selection, 'inside')
      % invert the selection, i.e., keep the voxels inside the box
      remove = ~remove;
    end

    if ismri
      if istrue(cfg.keepbrain)
        tmpcfg = [];
        tmpcfg.output = {'brain'};
        seg = ft_volumesegment(tmpcfg, mri);
        fprintf('keeping voxels in brain segmentation\n');
        % keep the tissue of the brain
        remove(seg.brain) = 0;
        clear seg
      end

      if istrue(cfg.feedback)
        tmpmri = keepfields(mri, {'anatomy', 'dim', 'transform', 'unit', 'coordsys'});
        tmpmri.remove = remove;
        tmpcfg = [];
        tmpcfg.funparameter = 'remove';
        ft_sourceplot(tmpcfg, tmpmri);
      end

      if isequal(cfg.smooth, 'no')
        fprintf('zero-filling %.0f%% of the volume\n', 100*mean(remove));
        mri.anatomy(remove) = 0;
      else
        tmp = mri.anatomy;
        tmp = (1 + 0.5.*randn(size(tmp))).*tmp; % add 50% noise to each voxel
        tmp = volumesmooth(tmp, cfg.smooth, 'anatomy');
        fprintf('smoothing %.0f%% of the volume\n', 100*mean(remove));
        mri.anatomy(remove) = tmp(remove);
      end

    elseif ismesh
      % determine all fields that might need to be defaced
      fn = setdiff(fieldnames(mri), ignorefields('deface'));
      dimord = cell(size(fn));
      for i=1:numel(fn)
        dimord{i} = getdimord(mri, fn{i});
      end
      % this applies to headshapes and meshes in general
      fprintf('keeping %d and removing %d vertices in the mesh\n', sum(remove==0), sum(remove==1));
      if isfield(mri, 'tri')
        [mri.pos, mri.tri] = remove_vertices(mri.pos, mri.tri, remove);
      elseif isfield(mri, 'tet')
        [mri.pos, mri.tet] = remove_vertices(mri.pos, mri.tet, remove);
      elseif isfield(mri, 'hex')
        [mri.pos, mri.hex] = remove_vertices(mri.pos, mri.hex, remove);
      else
        mri.pos = mri.pos(~remove,1:3);
      end
      for i=1:numel(fn)
        dimtok = tokenize(dimord{i}, '_');
        % do some sanity checks
        if any(strcmp(dimtok, '{pos}'))
          ft_error('not supported');
        end
        if numel(dimtok)>5
          ft_error('too many dimensions');
        end
        % remove the same positions from each matching dimension
        if numel(dimtok)>0 && strcmp(dimtok{1}, 'pos')
          mri.(fn{i}) = mri.(fn{i})(~remove,:,:,:,:);
        end
        if numel(dimtok)>1 && strcmp(dimtok{2}, 'pos')
          mri.(fn{i}) = mri.(fn{i})(:,~remove,:,:,:);
        end
        if numel(dimtok)>2 && strcmp(dimtok{3}, 'pos')
          mri.(fn{i}) = mri.(fn{i})(:,:,~remove,:,:);
        end
        if numel(dimtok)>3 && strcmp(dimtok{4}, 'pos')
          mri.(fn{i}) = mri.(fn{i})(:,:,:,~remove,:);
        end
        if numel(dimtok)>4 && strcmp(dimtok{5}, 'pos')
          mri.(fn{i}) = mri.(fn{i})(:,:,:,:,~remove);
        end
      end % for fn
      mri = removefields(mri, {'dim', 'transform'}); % these fields don't apply any more
    end % ismesh

  otherwise
    ft_error('unsupported method');
end % switch method

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble previous mri
ft_postamble provenance mri
ft_postamble history mri
ft_postamble savevar mri

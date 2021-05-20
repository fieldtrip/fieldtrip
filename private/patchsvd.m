function [sourcemodel] = patchsvd(cfg, sourcemodel)

% PATCHSVD computes a linear basis to span the leadfield for a defined patch
% of the source space. It is called by FT_PREPARE_LEADFIELD. This function 
% was originally written to do something like Limpiti et al.
% IEEE trans biomed eng 2006;53(9);1740-54, i.e. to create a linear basis 
% to span the leadfield for a patch of cortex, based on an SVD. It now also
% implements the procedure to compute a (spatial basis) for a ROI's
% leadfield, e.g. as per Backus et al. DOI:10.1016/j.cub.2015.12.048.
%
% Supported cfg options are:
%   cfg.patchsvd = 'yes', or a scalar. The scalar value is to support old
%                  behavior, in which case it is treated as a distance to 
%                  define the inclusion of dipoles to define the patch
%   cfg.patchsvdnum = scalar, integer number or percentage, defining the
%                  number of spatial components per patch, or the total
%                  amount of 'spatial variance' explained by the the
%                  patch' basis. Default is 5.
%   cfg.atlas    = a specification of an atlas to be used for the
%                  definition of the patches
%           
%   cfg.parcellation  = string, name of the atlas field that is used for the
%                       parcel labels. (default = [])
%   cfg.parcel        = string, or cell-array of strings, specifying for which
%                       parcels to return the output. (default = 'all')
%
% See also FT_VIRTUALCHANNEL

% Copyright (c) 2006-2021, Jan-Mathijs Schoffelen & Robert Oostenveld, DCCN
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

% set the defaults
cfg.patchsvd    = ft_getopt(cfg, 'patchsvd',    3); % this is defined as a distance
cfg.patchsvdnum = ft_getopt(cfg, 'patchsvdnum', 5); % this is defined as a number of spatial components to be kept
cfg.feedback    = ft_getopt(cfg, 'feedback',    'textbar');

hasatlas = isfield(cfg, 'atlas');
if hasatlas && (isnumeric(cfg.patchsvd) || isequal(cfg.patchsvd, 'all'))
  ft_error('the specification of an atlas is incompatible with a numeric patchsvd in the cfg, or with cfg.patchsvd = ''all''');
end

if hasatlas
  cfg.parcellation = ft_getopt(cfg, 'parcellation');
  cfg.parcel       = ft_getopt(cfg, 'parcel',   'all');

  % ensure it is a parcellation, not a segmentation
  parcellation = ft_checkdata(cfg.atlas, 'datatype', 'parcellation', 'parcellationstyle', 'indexed');
  
  % ensure that the source and the parcellation are anatomically consistent
  if ~isalmostequal(sourcemodel.pos, parcellation.pos, 'abstol', 1000000*eps)
    ft_error('the source positions are not consistent with the parcellation, please use FT_SOURCEINTERPOLATE');
  end
  
  if isequal(cfg.parcel, 'all')
    cfg.parcel = parcellation.(sprintf('%slabel',cfg.parcellation));
  end
  
  indx       = match_str(parcellation.(sprintf('%slabel',cfg.parcellation)), cfg.parcel);
  cfg.parcel = parcellation.(sprintf('%slabel',cfg.parcellation))(indx);
  
  n     = zeros(numel(indx),1);
  lfall = cell(numel(indx),1);
  pos   = zeros(numel(indx),3);
  ft_progress('init', cfg.feedback, 'computing patchsvd');
  for k = 1:numel(indx)
    ft_progress(indx(k)/numel(indx), 'computing patchsvd for parcel %s\n', cfg.parcel{k});
    
    sel      = find(parcellation.(cfg.parcellation)==indx(k));
    pos(k,:) = mean(sourcemodel.pos(sel,:),1);
    
    % concatenate the leadfield of all dipoles that are inside the ROI into one matrix
    lfr     = sourcemodel.leadfield(sel(:));
    lfr     = cat(2, lfr{:});
    % svd of leadfields of dipoles inside the ROI
    [U,S,V] = svd(lfr, 'econ');
    
    if isequal(cfg.patchsvdnum, 'all')
      n(k) = size(S,1);
    elseif cfg.patchsvdnum < 1
      % Limpiti et al 2006 formula 12
      s    = diag(S).^2;
      s    = cumsum(s)./sum(s);
      n(k) = find(s - cfg.patchsvdnum > 0, 1);
    else
      n(k) = min(numel(sel), cfg.patchsvdnum);
    end
    %lfall{k} = U(:,1:n(k)); %*S(1:n(dipindx),1:n(dipindx)); %klopt dit?
    %coeff{k} = V(1:n(k), :);
    lfall{k} = lfr * V(:,1:n(k));
  end
  ft_progress('close');
  
  sourcemodel     = removefields(sourcemodel, {'pos' 'dim' 'tri' 'inside'});
  sourcemodel.pos = pos; 
  sourcemodel.inside = true(size(pos,1),1);
  
elseif isnumeric(cfg.patchsvd)
  % this is based on the old (probably never used) code 
  Ndipoles = size(sourcemodel.pos,1);
  Ninside  = sum(sourcemodel.inside);
  inside   = find(sourcemodel.inside);
  
  lfall    = cell(Ndipoles,1);
  n        = zeros(Ndipoles,1);
  
  ft_progress('init', cfg.feedback, 'computing patchsvd');
  for dipindx=1:Ninside
    % renumber the loop-index variable to make it easier to print the progress bar
    i = inside(dipindx);
    
    % compute the distance from this dipole to each other dipole
    dist = sqrt(sum((sourcemodel.pos-repmat(sourcemodel.pos(i,:), [Ndipoles 1])).^2, 2));
    
    % define the region of interest around this dipole
    sel  = dist<=cfg.patchsvd & sourcemodel.inside;
    Nsel = sum(sel);
    
    ft_progress(dipindx/Ninside, 'computing patchsvd %d/%d, Nsel=%d\n', dipindx, Ninside, Nsel);
    % concatenate the leadfield of all dipoles that are inside the ROI into one matrix
    lfr     = sourcemodel.leadfield(sel(:)');
    lfr     = cat(2,lfr{:});
    % svd of leadfields of dipoles inside the ROI
    [U,S,V] = svd(lfr, 'econ');
    
    if cfg.patchsvdnum < 1
      % Limpiti et al 2006 formula 12
      s          = diag(S).^2;
      s          = cumsum(s)./sum(s);
      n(dipindx) = find(s - cfg.patchsvdnum > 0, 1);
    else
      n(dipindx) = cfg.patchsvdnum;
    end
    lfall{i} = lfr*V(:,1:n(dipindx));% U(:,1:n(dipindx)); %*S(1:n(dipindx),1:n(dipindx)); %klopt dit?
    %nghbr{i} = sel;
    %coeff{i} = V(1:n(dipindx), :);
  end
  ft_progress('close');
  
elseif isequal(cfg.patchsvd, 'all')
  lfr     = sourcemodel.leadfield(sourcemodel.inside(:)');
  lfr     = cat(2, lfr{:});
  [U,S,V] = svd(lfr, 'econ');

  if cfg.patchsvdnum < 1
    % Limpiti et al 2006 formula 12
    s = diag(S).^2;
    s = cumsum(s)./sum(s);
    n = find(s - cfg.patchsvdnum > 0, 1);
  else
    n = cfg.patchsvdnum;
  end
  lfall{1} = lfr*V(:,1:n); % U(:,1:n); %*S(1:n(dipindx),1:n(dipindx)); %klopt dit?
  %nghbr{1} = sourcemodel.inside;
  %coeff{1} = V(1:n, :);

  %---change output
  sourcemodel.pos    = mean(sourcemodel.pos(sourcemodel.inside,:),1);
  sourcemodel.inside = true(1);
  sourcemodel.dim    = [1 1 1];
else
  %do nothing
end

% update the leadfields
sourcemodel.leadfield = lfall;

function bnd = ft_iso2surf(cfg, mri)

% FT_ISO2SURF is an alternative to ft_prepare_mesh for generating surfaces 
% from tissue-labeled MRIs. Requires the ISO2MESH matlab toolbox.
% cfg.tissue          = cell-array with tissue types (for segmentation
%                       data) or numeric vector with integer values (for 
%                       tissue-labeled MRI)
% cfg.numvertices     = numeric vector, should have same number of elements
%                       as cfg.tissue
% cfg.ascending       = indicates whether tissue labeled MRI values are 
%                       ascending from outer-inner layers. Only needed for
%                       tissue-labeled MRIs. ['yes']/'no'
% mri                 = tissue-labeled MRI or segmentation data structure
%
% Copyright (C) 2013, Daniel D.E. Wong, Sarang S. Dalal
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
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%

if(~exist('vol2surf','file'))
    error('The iso2mesh library is required in the MATLAB path. It can be downloaded from http://iso2mesh.sourceforge.net');
end

%mri = ft_checkdata(mri, 'datatype', {'volume', 'segmentation'});
%mri = ft_convert_units(mri);

% basedonmri        = ft_datatype(mri, 'volume');
basedonseg        = ft_datatype(mri, 'segmentation');

if ~basedonseg
    assert(isnumeric(cfg.tissue), 'cfg.tissue should contain numeric tissue label values for tissue-labeled MRI data.');
    
    if isfield(cfg,'ascending') && strcmp(cfg.ascending,'no')
        mri.anatomy = -mri.anatomy;
        cfg.tissue = -cfg.tissue;
    end
end

for ii = 1:length(cfg.tissue)
    if ~basedonseg
        Ymesh = Y>=cfg.tissue(ii);
    else
        Ymesh = mri.(cfg.tissue{ii});
    end
    
    opt.radbound=3;    % set the target surface mesh element bounding sphere be <3 pixels in radius
    opt.maxnode=cfg.numvertices(ii);
    opt.maxsurf=1;
    [node,face]=v2s(Ymesh,1,opt,'cgalsurf');
    face = face(:,1:3);
    %face = meshreorient(node,face);
    
    % Convert voxels indexes as expected by SPM - this needs testing
%     d = sign(diag(mri.transform))';
%     node = node.*repmat(d(1:3),size(node,1),1);
%     for jj = 1:3; node(:,jj) = mod(node(:,jj),mri.dim(jj)); end; %.*abs(d(jj))); end;
    node = [node ones(size(node,1),1)] * mri.transform';
    node = node(:,1:3);
    
    bnd(ii).pnt = node;
    bnd(ii).tri = face;
    bnd(ii).unit = mri.unit;
    region_size(ii) = sum(find(Ymesh));
end

% Order outside in (trying to be smart here)
[~,order] = sort(region_size,'descend');
bnd = bnd(order);

bnd = decouplesurf(bnd);

for ii = 1:length(bnd)
    [bnd(ii).pnt, bnd(ii).tri] = surfreorient(bnd(ii).pnt, bnd(ii).tri);
    bnd(ii).tri = bnd(ii).tri(:,[3 2 1]);
end


%%
% Take care of intersecting surfaces.
%%
function bnd = decouplesurf(bnd)

for ii = 1:length(bnd)-1
    % Despite what the instructions for surfboolean says, surfaces should
    % be ordered from inside-out!!
    [newnode, newelem] = surfboolean(bnd(ii+1).pnt,bnd(ii+1).tri,'decouple',bnd(ii).pnt,bnd(ii).tri);
    bnd(ii+1).tri = newelem(newelem(:,4)==2,1:3) - size(bnd(ii+1).pnt,1);
    bnd(ii+1).pnt = newnode(newnode(:,4)==2,1:3);
end
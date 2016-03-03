function [headmodel, cfg] = ft_prepare_bemmodel(cfg, mri)

% FT_PREPARE_BEMMODEL  is deprecated, please use FT_PREPARE_HEADMODEL and
% FT_PREPARE_MESH
%
% See also FT_PREPARE_HEADMODEL

% Copyright (C) 2005-2009, Robert Oostenveld
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

warning('FT_PREPARE_BEMMODEL is deprecated, please use FT_PREPARE_HEADMODEL with cfg.method = ''dipoli/openmeeg/bemcp'' instead.')

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble provenance mri
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

cfg = ft_checkconfig(cfg, 'renamed', {'hdmfile', 'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'vol',     'headmodel'});

% set the defaults
if ~isfield(cfg, 'tissue'),         cfg.tissue = [8 12 14];                  end
if ~isfield(cfg, 'numvertices'),    cfg.numvertices = [1 2 3] * 500;         end
if ~isfield(cfg, 'isolatedsource'), cfg.isolatedsource = [];                 end
if ~isfield(cfg, 'method'),         cfg.method = 'dipoli';                   end % dipoli, openmeeg, bemcp

% start with an empty volume conductor
try
  hdm = ft_fetch_vol(cfg);
  headmodel.bnd = hdm.bnd;
  if isfield(hdm, 'cond')
    % also copy the conductivities
    headmodel.cond = hdm.cond;
  end
catch
  headmodel = [];
  geom = fixpos(mri);
  % copy the boundaries from the geometry into the volume conduction model
  headmodel.bnd = geom.bnd;
end

% determine the number of compartments
Ncompartment = numel(headmodel.bnd);

% assign the conductivity
if ~isfield(headmodel,'cond')
  if ~isfield(cfg, 'conductivity')
    if isfield(mri, 'cond')
      headmodel.cond = mri.cond;
    elseif isfield(mri, 'c')
      headmodel.cond = mri.c;
    else
      fprintf('warning: using default values for the conductivity')
      headmodel.cond = [1 1/80 1] * 0.33;
    end
  else
    if ~isempty(cfg.conductivity)
      headmodel.cond = cfg.conductivity;
    elseif isempty(cfg.conductivity) && Ncompartment==3
      fprintf('warning: using default values for the conductivity')
      headmodel.cond = [1 1/80 1] * 0.33;
    else
      fprintf('warning: using 1 for all conductivities')
      headmodel.cond = ones(1,Ncompartment);
    end
  end
end

if ~isfield(headmodel, 'bnd')
  % construct the geometry of the BEM boundaries
  if nargin==1
    headmodel.bnd = ft_prepare_mesh(cfg);
  else
    headmodel.bnd = ft_prepare_mesh(cfg, mri);
  end
end

headmodel.source = find_innermost_boundary(headmodel.bnd);
headmodel.skin_surface   = find_outermost_boundary(headmodel.bnd);
fprintf('determining source compartment (%d)\n', headmodel.source);
fprintf('determining skin compartment (%d)\n',   headmodel.skin_surface);

if ~isempty(cfg.isolatedsource)
  isolatedsource = istrue(cfg.isolatedsource);
else
  isolatedsource = false;
end

if isempty(cfg.isolatedsource) && Ncompartment>1 && strcmp(cfg.method, 'dipoli')
  % the isolated source compartment is by default the most inner one
  isolatedsource = true;
elseif isempty(cfg.isolatedsource) && Ncompartment==1
  % the isolated source interface should be contained within at least one other interface
  isolatedsource = false;
elseif ~islogical(isolatedsource)
  error('cfg.isolatedsource should be true or false');
end

if cfg.isolatedsource
  fprintf('using compartment %d for the isolated source approach\n', headmodel.source);
else
  fprintf('not using the isolated source approach\n');
end

if strcmp(cfg.method, 'dipoli')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % this uses an implementation that was contributed by Thom Oostendorp
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ft_hastoolbox('dipoli', 1);

  % use the dipoli wrapper function
  headmodel = dipoli(headmodel, isolatedsource);
  headmodel.type = 'dipoli';

elseif strcmp(cfg.method, 'bemcp')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % this uses an implementation that was contributed by Christophe Philips
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ft_hastoolbox('bemcp', 1);

  % do some sanity checks
  if length(headmodel.bnd)~=3
    error('this only works for three surfaces');
  end

  numboundaries = length(headmodel.bnd);
  % determine the nesting of the compartments
  nesting = zeros(numboundaries);
  for i=1:numboundaries
    for j=1:numboundaries
      if i~=j
        % determine for a single vertex on each surface if it is inside or outside the other surfaces
        curpos = headmodel.bnd(i).pos(1,:); % any point on the boundary is ok
        curpnt = headmodel.bnd(j).pos;
        curtri = headmodel.bnd(j).tri;
        nesting(i,j) = bounding_mesh(curpos, curpnt, curtri);
      end
    end
  end

  if sum(nesting(:))~=(numboundaries*(numboundaries-1)/2)
    error('the compartment nesting cannot be determined');
  end

  % for a three compartment model, the nesting matrix should look like
  %    0 1 1     the first is nested inside the 2nd and 3rd, i.e. the inner skull
  %    0 0 1     the second is nested inside the 3rd, i.e. the outer skull
  %    0 0 0     the third is the most outside, i.e. the skin
  [dum, order] = sort(-sum(nesting,2));

  fprintf('reordering the boundaries to: ');
  fprintf('%d ', order);
  fprintf('\n');

  % update the order of the compartments
  headmodel.bnd    = headmodel.bnd(order);
  headmodel.cond   = headmodel.cond(order);
  headmodel.skin_surface   = numboundaries;
  headmodel.source = 1;

  % Build Triangle 4th point
  headmodel = triangle4pt(headmodel);

  % 2. BEM model estimation, only for the scalp surface

  defl =[ 0 0 1/size(headmodel.bnd(headmodel.skin_surface).pos,1)];
  % ensure deflation for skin surface, i.e. average reference over skin

  % NOTE:
  % Calculation proceeds by estimating each submatrix C_ij and combine them.
  % There are 2 options:
  % - calculating the matrices once, as it takes some time, keep them in
  %   memory and use them the 2-3 times they're needed.
  % - calculating the matrices every time they're needed, i.e. 2-3 times
  % The former option is faster but requires more memory space as up to *8*
  % square matrices of size C_ij have to be kept in memory at once.
  % The latter option requires less memory, but would take much more time to
  % estimate.
  % This faster but memory hungry solution is implemented here.

  % Deal first with surface 1 and 2 (inner and outer skull
  %--------------------------------

  % NOTE:
  % C11st/C22st/C33st are simply the matrix C11/C22/C33 minus the identity
  % matrix, i.e. C11st = C11-eye(N)

  weight = (headmodel.cond(1)-headmodel.cond(2))/((headmodel.cond(1)+headmodel.cond(2))*2*pi);
  C11st  = bem_Cii_lin(headmodel.bnd(1).tri,headmodel.bnd(1).pos, weight,defl(1),headmodel.bnd(1).pnt4);
  weight = (headmodel.cond(1)-headmodel.cond(2))/((headmodel.cond(2)+headmodel.cond(3))*2*pi);
  C21    = bem_Cij_lin(headmodel.bnd(2).pos,headmodel.bnd(1).pos,headmodel.bnd(1).tri, weight,defl(1));
  tmp1   = C21/C11st;

  weight = (headmodel.cond(2)-headmodel.cond(3))/((headmodel.cond(1)+headmodel.cond(2))*2*pi);
  C12    = bem_Cij_lin(headmodel.bnd(1).pos,headmodel.bnd(2).pos,headmodel.bnd(2).tri, weight,defl(2));
  weight = (headmodel.cond(2)-headmodel.cond(3))/((headmodel.cond(2)+headmodel.cond(3))*2*pi);
  C22st  = bem_Cii_lin(headmodel.bnd(2).tri,headmodel.bnd(2).pos, weight,defl(2),headmodel.bnd(2).pnt4);
  tmp2   = C12/C22st;

  % Try to spare some memory:
  tmp10 = - tmp2 * C21 + C11st;
  clear C21 C11st
  tmp11 = - tmp1 * C12 + C22st;
  clear C12 C22st

  % Combine with the effect of surface 3 (scalp) on the first 2
  %------------------------------------------------------------
  weight = (headmodel.cond(1)-headmodel.cond(2))/(headmodel.cond(3)*2*pi);
  C31    = bem_Cij_lin(headmodel.bnd(3).pos,headmodel.bnd(1).pos,headmodel.bnd(1).tri, weight,defl(1));
  %   tmp4   = C31/(- tmp2 * C21 + C11st );
  %   clear C31 C21 C11st
  tmp4 = C31/tmp10;
  clear C31 tmp10

  weight = (headmodel.cond(2)-headmodel.cond(3))/(headmodel.cond(3)*2*pi);
  C32    = bem_Cij_lin(headmodel.bnd(3).pos,headmodel.bnd(2).pos,headmodel.bnd(2).tri, weight,defl(2));
  %   tmp3   = C32/(- tmp1 * C12 + C22st );
  %   clear  C12 C22st C32
  tmp3 = C32/tmp11;
  clear C32 tmp11

  tmp5 = tmp3*tmp1-tmp4;
  tmp6 = tmp4*tmp2-tmp3;
  clear tmp1 tmp2 tmp3 tmp4

  % Finally include effect of surface 3 on the others
  %--------------------------------------------------
  % As the gama1 intermediate matrix is built as the sum of 3 matrices, I can
  % spare some memory by building them one at a time, and summing directly
  weight = headmodel.cond(3)/((headmodel.cond(1)+headmodel.cond(2))*2*pi);
  Ci3    = bem_Cij_lin(headmodel.bnd(1).pos,headmodel.bnd(3).pos,headmodel.bnd(3).tri, weight,defl(3));
  gama1  = - tmp5*Ci3; % gama1 = - tmp5*C13;

  weight = headmodel.cond(3)/((headmodel.cond(2)+headmodel.cond(3))*2*pi);
  Ci3    = bem_Cij_lin(headmodel.bnd(2).pos,headmodel.bnd(3).pos,headmodel.bnd(3).tri, weight,defl(3));
  gama1  = gama1 - tmp6*Ci3; % gama1 = - tmp5*C13 - tmp6*C23;

  weight = 1/(2*pi);
  Ci3    = bem_Cii_lin(headmodel.bnd(3).tri,headmodel.bnd(3).pos, weight,defl(3),headmodel.bnd(3).pnt4);
  gama1  = gama1 - Ci3; % gama1 = - tmp5*C13 - tmp6*C23 - C33st;
  clear Ci3

  % Build system matrix
  %--------------------
  i_gama1 = inv(gama1);
  headmodel.mat = [i_gama1*tmp5 i_gama1*tmp6 i_gama1];
  headmodel.type = 'bemcp';

elseif strcmp(cfg.method, 'openmeeg')
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % this uses an implementation that was contributed by INRIA Odyssee Team
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ft_hastoolbox('openmeeg', 1);

  if size(headmodel.bnd(1).pos,1)>10000
    error('OpenMEEG does not manage meshes with more than 10000 vertices (use reducepatch)')
  else
    % use the openmeeg wrapper function
    headmodel = openmeeg(headmodel,cfg.isolatedsource);
    headmodel.type = 'openmeeg';
  end

else
  error('unsupported method');
end % which method

% ensure that the geometrical units are specified
headmodel = ft_convert_units(headmodel);

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
ft_postamble provenance headmodel
ft_postamble history    headmodel

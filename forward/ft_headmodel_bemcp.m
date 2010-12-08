function vol = ft_headmodel_bemcp(geom, varargin)

% FT_HEADMODEL_BEMCP creates a volume conduction model of the head
% using the boundary element method (BEM) for EEG. This function
% takes as input the triangulated surfaces that describe the boundaries
% and returns as output a volume conduction model which can be used
% to compute leadfields.
%
% The implementation of this function is based on Christophe Phillips'
% MATLAB code, hence the name "bemcp".
%
% Use as
%   vol = ft_headmodel_bem_cp(geom, ...)
%
% See also FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD

ft_hastoolbox('bemcp', 1);

% get the optional arguments
hdmfile         = keyval('hdmfile', varargin);
conductivity    = keyval('conductivity', varargin);

% start with an empty volume conductor
vol = [];

if ~isempty(hdmfile)
  hdm = ft_read_vol(hdmfile);
  % copy the boundary of the head model file into the volume conduction model
  vol.bnd = hdm.bnd;
  if isfield(hdm, 'cond')
    % also copy the conductivities
    vol.cond = hdm.cond;
  end
else
  % copy the boundaries from the geometry into the volume conduction model
  vol.bnd = geom.bnd;
end

% determine the number of compartments
numboundaries = length(vol.bnd);

if ~isfield(vol, 'cond')
  % assign the conductivity of each compartment
  vol.cond = conductivity;
end

% determine the nesting of the compartments
nesting = zeros(numboundaries);
for i=1:numboundaries
  for j=1:numboundaries
    if i~=j
      % determine for a single vertex on each surface if it is inside or outside the other surfaces
      curpos = vol.bnd(i).pnt(1,:); % any point on the boundary is ok
      curpnt = vol.bnd(j).pnt;
      curtri = vol.bnd(j).tri;
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
[~, order] = sort(-sum(nesting,2));

fprintf('reordering the boundaries to: ');
fprintf('%d ', order);
fprintf('\n');

% update the order of the compartments
vol.bnd    = vol.bnd(order);
vol.cond   = vol.cond(order);
vol.skin_surface   = numboundaries;
vol.source = 1;

% do some sanity checks
if length(vol.bnd)~=3
  error('this only works for three surfaces');
end
if vol.skin_surface~=3
  error('the third surface should be the skin');
end
if vol.source~=1
  error('the first surface should be the inside of the skull');
end

% Build Triangle 4th point
vol = triangle4pt(vol);

% 2. BEM model estimation, only for the scalp surface

defl =[ 0 0 1/size(vol.bnd(vol.skin_surface).pnt,1)];
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

weight = (vol.cond(1)-vol.cond(2))/((vol.cond(1)+vol.cond(2))*2*pi);
C11st  = bem_Cii_lin(vol.bnd(1).tri,vol.bnd(1).pnt, weight,defl(1),vol.bnd(1).pnt4);
weight = (vol.cond(1)-vol.cond(2))/((vol.cond(2)+vol.cond(3))*2*pi);
C21    = bem_Cij_lin(vol.bnd(2).pnt,vol.bnd(1).pnt,vol.bnd(1).tri, weight,defl(1));
tmp1   = C21/C11st;

weight = (vol.cond(2)-vol.cond(3))/((vol.cond(1)+vol.cond(2))*2*pi);
C12    = bem_Cij_lin(vol.bnd(1).pnt,vol.bnd(2).pnt,vol.bnd(2).tri, weight,defl(2));
weight = (vol.cond(2)-vol.cond(3))/((vol.cond(2)+vol.cond(3))*2*pi);
C22st  = bem_Cii_lin(vol.bnd(2).tri,vol.bnd(2).pnt, weight,defl(2),vol.bnd(2).pnt4);
tmp2   = C12/C22st;

% Try to spare some memory:
tmp10 = - tmp2 * C21 + C11st;
clear C21 C11st
tmp11 = - tmp1 * C12 + C22st;
clear C12 C22st

% Combine with the effect of surface 3 (scalp) on the first 2
%------------------------------------------------------------
weight = (vol.cond(1)-vol.cond(2))/(vol.cond(3)*2*pi);
C31    = bem_Cij_lin(vol.bnd(3).pnt,vol.bnd(1).pnt,vol.bnd(1).tri, weight,defl(1));
%   tmp4   = C31/(- tmp2 * C21 + C11st );
%   clear C31 C21 C11st
tmp4 = C31/tmp10;
clear C31 tmp10

weight = (vol.cond(2)-vol.cond(3))/(vol.cond(3)*2*pi);
C32    = bem_Cij_lin(vol.bnd(3).pnt,vol.bnd(2).pnt,vol.bnd(2).tri, weight,defl(2));
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
weight = vol.cond(3)/((vol.cond(1)+vol.cond(2))*2*pi);
Ci3    = bem_Cij_lin(vol.bnd(1).pnt,vol.bnd(3).pnt,vol.bnd(3).tri, weight,defl(3));
gama1  = - tmp5*Ci3; % gama1 = - tmp5*C13;

weight = vol.cond(3)/((vol.cond(2)+vol.cond(3))*2*pi);
Ci3    = bem_Cij_lin(vol.bnd(2).pnt,vol.bnd(3).pnt,vol.bnd(3).tri, weight,defl(3));
gama1  = gama1 - tmp6*Ci3; % gama1 = - tmp5*C13 - tmp6*C23;

weight = 1/(2*pi);
Ci3    = bem_Cii_lin(vol.bnd(3).tri,vol.bnd(3).pnt, weight,defl(3),vol.bnd(3).pnt4);
gama1  = gama1 - Ci3; % gama1 = - tmp5*C13 - tmp6*C23 - C33st;
clear Ci3

% Build system matrix
%--------------------
i_gama1 = inv(gama1);
vol.mat = [i_gama1*tmp5 i_gama1*tmp6 i_gama1];

% remember the type
vol.type = 'bemcp';

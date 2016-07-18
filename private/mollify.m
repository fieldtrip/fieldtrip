function [grid] = mollify(cfg, grid)

% This function does something
%

% Copyright (c) 2006, Jan-Mathijs Schoffelen & Robert Oostenveld, F.C. Donders Centre
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
if ~isfield(cfg, 'mollify'),     cfg.mollify = 1;           end  % fwhm in pos-units
if ~isfield(cfg, 'sphereradius'),cfg.sphereradius = 2.*cfg.mollify; end % truncate gaussian 
if ~isfield(cfg, 'feedback'),    cfg.feedback = 'textbar';     end

hasnrm   = isfield(grid, 'normals');

Ndipoles = size(grid.pos,1);
Ninside  = length(grid.inside);
Nchans   = size(grid.leadfield{grid.inside(1)}, 1);
Ncomp    = size(grid.leadfield{grid.inside(1)}, 2);

if isempty(cfg.sphereradius)
  error('cfg.sphereradius should be specified');
end
% the distance only has to be computed to voxels inside the brain
pos  = grid.pos(grid.inside,:);
npos = size(pos,1);
% compute the distance between voxels and each surface point in a reasonable amount of time, so don't compute everything with everything, but only take those point-pairs into account at a distance of < sphereradius, because the gaussian will be truncated anyhow
distmat      = pntdist(pos, cfg.sphereradius);
[indx, indy] = find(distmat > 0 & distmat <= cfg.sphereradius);
val          = distmat(find(distmat > 0 & distmat <= cfg.sphereradius));
distmat      = sparse(indx, indy, val, Ninside, Ninside);

ldfall   = zeros(Nchans*Ncomp,Ninside);
%put leadfields in one matrix and concatenate the columns
for k=1:Ninside
  ldfall(:,k) = grid.leadfield{grid.inside(k)}(:);
end

if hasnrm,
  nrmall = grid.normals(grid.inside,:)';
  nrmnew = zeros(size(nrmall));
end

ldfnew   = cell(1,Ndipoles);
sigma    = cfg.mollify./(2*sqrt(2*log(2))); %Weisstein, Eric W. "Gaussian Function." From MathWorld--A Wolfram Web Resource. http://mathworld.wolfram.com/GaussianFunction.html 

progress('init', cfg.feedback, 'computing mollification');
for k=1:Ninside
  progress(k/Ninside, 'computing mollification %d/%d\n', k, Ninside);
  % compute the squared distance from this dipole to each other dipole
  distsq = full(distmat(:,k).^2);
  % gaussianize the kernel
  kernel = exp(-distsq./(2.*(sigma^2)));%CHECK THIS
  % put everything outside the sphereradius to zero, except the point itself
  kernel(find(distsq==0)) = 0;
  kernel(k) = 1; 
  % normalize kernel
  kernel = kernel./sum(kernel);
  % compute mollified leadfield
  dum = reshape(ldfall*kernel, [Nchans Ncomp]);
  if ~strcmp(cfg.reducerank, 'no'),
    [u,s,v] = svd(dum);
    r       = diag(s);
    s(:)    = 0;
    s(1,1)  = r(1);
    s(2,2)  = r(2);
    dum     = u*s*v';
  end
  ldfnew{grid.inside(k)} = dum;
  if hasnrm, nrmnew(:, grid.inside(k)) = nrmall * kernel; end
end
progress('close');
  
% update the leadfields and the normals
grid.leadfield = ldfnew;
if hasnrm, grid.normals = nrmnew; end

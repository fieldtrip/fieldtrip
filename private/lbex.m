function [grid] = lbex(cfg, grid)

% This function will add the field "subspace" to the grid definition.
%
% The subspace projection is based on the LBEX (local basis expansion) 
% method.

% set the defaults
if ~isfield(cfg, 'lbex'),       cfg.lbex = 3;              end
if ~isfield(cfg, 'lbexeigtol'), cfg.lbexeigtol = 1000*eps; end
if ~isfield(cfg, 'feedback'),   cfg.feedback = 'text';     end

Ndipoles = size(grid.pos,1);
Ninside  = length(grid.inside);

% concatenate the leadfield of all dipoles that are inside the brain into one large matrix
sel = grid.inside;
lfa = cell2mat(grid.leadfield(sel(:)'));
% covariance of all leadfields
Ca = lfa * lfa';

progress('init', cfg.feedback, 'computing lbex');
for dipindx=1:Ninside
  % renumber the loop-index variable to make it easier to print the progress bar
  i = grid.inside(dipindx);
  
  % compute the distance from this dipole to each other dipole
  dist = sqrt(sum((grid.pos-repmat(grid.pos(i,:), [Ndipoles 1])).^2, 2));

  % define the region of interest around this dipole
  sel  = find(dist<=cfg.lbex);
  sel  = intersect(sel, grid.inside);
  Nsel = length(sel);

  progress(dipindx/Ninside, 'computing lbex %d/%d, Nsel=%d\n', dipindx, Ninside, Nsel);

  % concatenate the leadfield of all dipoles that are inside the ROI into one matrix
  lfr = cell2mat(grid.leadfield(sel(:)'));
  % covariance of leadfields of dipoles inside the ROI
  Cr = lfr * lfr';
  
  % The eigenvalue problem is to determine the nontrivial solutions of the equation
  %   A*x = l*x
  % The generalized eigenvalue decomposition solves
  %   A*x = l*B*x
  % If B is non-singular, the problem could be solved by reducing it into a standard eigenvalue problem 
  %   inv(B)*A*x = l*x 
  % See http://www.mathworks.com/access/helpdesk/help/techdoc/ref/eig.html

  % compute eigenspace decomposition, THIS IS NUMERICALLY UNSTABLE
  [v, d] = eig(Cr, Ca);
  % compute eigenspace decomposition, THIS ALSO DOES NOT SOLVE THE PROBLEM
  % [v, d] = eig(pinv(Ca, cfg.lbexeigtol)*Cr);

  % select the eigenvectors with a non-zero eigenvalue
  dd  = diag(d);
  dd  = dd./max(dd);
  sel = dd>cfg.lbexeigtol;
  
  % remember the subspace projection matrix
  grid.subspace{grid.inside(dipindx)} = v(:, sel)';
end
progress('close');

% fill the positions outside the brain with NaNs
for dipindx=grid.outside(:)'
  grid.subspace{dipindx} = nan;
end


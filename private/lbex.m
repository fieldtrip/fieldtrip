function [sourcemodel] = lbex(cfg, sourcemodel)

% This function will add the field "subspace" to the sourcemodel definition.
%
% The subspace projection is based on the LBEX (local basis expansion)
% method.

% set the defaults
cfg.lbex       = ft_getopt(cfg, 'lbex',       3); % this is a distance, in units of sourcemodel.pos
cfg.lbexeigtol = ft_getopt(cfg, 'lbexeigtol', 1000*eps);
cfg.feedback   = ft_getopt(cfg, 'feedback',   'text');
cfg.keep       = ft_getopt(cfg, 'keep',       'all');

if isequal(cfg.keep, 'all')
  cfg.keep = sourcemodel.inside;
else
  if ~islogical(cfg.keep)
    keep = false(size(sourcemodel.pos,1),1);
    keep(cfg.keep) = true;
    cfg.keep = keep;
  end
end
assert(isequal(numel(cfg.keep), numel(sourcemodel.inside)));

cfg.keep = cfg.keep(:) & sourcemodel.inside(:);

Ndipoles = size(sourcemodel.pos,1);
Ninside  = sum(cfg.keep);
inside   = find(cfg.keep);

% concatenate the leadfield of all dipoles that are inside the brain into one large matrix
lfa = cat(2, sourcemodel.leadfield{:});

% do the computations on the svd basis to avoid numerical issues
[U,S,V] = svd(lfa, 'econ');
diagS   = diag(S);
Tol     = 1e-12;
sel     = find(diagS>Tol.*diagS(1));
P       = diag(1./sqrt(diag(S(sel,sel))))*U(:,sel)'; % prewhitening matrix
lfa     = P*lfa;

% covariance of all leadfields
Ca = lfa * lfa';

sourcemodel.subspace = cell(1,size(sourcemodel.pos,1));
ft_progress('init', cfg.feedback, 'computing lbex');
for dipindx=1:Ninside
  % renumber the loop-index variable to make it easier to print the progress bar
  i = inside(dipindx);

  % compute the distance from this dipole to each other dipole
  dist = sqrt(sum((sourcemodel.pos-repmat(sourcemodel.pos(i,:), [Ndipoles 1])).^2, 2));
    
  % define the region of interest around this dipole
  sel  = dist<=cfg.lbex & sourcemodel.inside;
  Nsel = sum(sel);

  % concatenate the leadfield of all dipoles that are inside the ROI into one matrix
  lfr = P*cat(2,sourcemodel.leadfield{sel(:)'});
  
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
  Nsel2 = sum(sel);
  
  ft_progress(dipindx/Ninside, 'computing lbex %d/%d, number of dipoles in ROI=%d, subspace dimension=%d\n', dipindx, Ninside, Nsel, Nsel2);

  % remember the subspace projection matrix
  sourcemodel.subspace{inside(dipindx)} = flip(v(:, sel)',1)*P;
end
ft_progress('close');

if ~isequal(cfg.keep(:),sourcemodel.inside(:))
  sourcemodel.inside = cfg.keep;
  sourcemodel.leadfield(~cfg.keep) = {[]};
  sourcemodel.subspace(~cfg.keep)  = {[]};
end

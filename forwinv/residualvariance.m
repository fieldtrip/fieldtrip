function [dipout] = residualvariance(dip, grad, vol, dat, varargin);

% RESIDUALVARIANCE scan with a single dipole and computes the RV
% at each grid location
% 
% Use as 
%   [dipout] = residualvariance(dip, grad, vol, dat, ...)

% Copyright (C) 2004-2006, Robert Oostenveld
% 
% $Log: residualvariance.m,v $
% Revision 1.5  2006/05/10 08:18:21  roboos
% swiched to using keyval() function for getting optional arguments instead of using eval()
%
% Revision 1.4  2005/10/25 08:55:24  roboos
% added some feedback for subspace projection
% changed indentation and some whitespace
%
% Revision 1.3  2004/10/27 16:16:37  roboos
% renamed grid.lbex matrix for subspace projection into grid.subspace (in correspondence with precompute_leadfield)
% transposed the subspace projection matrix
% renamed all occurences of "lbex" (as part of variable names) into "subspace"
%
% Revision 1.2  2004/10/25 16:22:21  roboos
% fixed parsing of optional arguments (key,value-pairs)
%
% Revision 1.1  2004/09/28 14:32:29  roboos
% initial version, implements LBEX
%

% get the optional settings, or use default value
feedback      = keyval('feedback',      varargin); if isempty(feedback),      feedback = 'text';            end

% ensure that these are row-vectors
dip.inside = dip.inside(:)';
dip.outside = dip.outside(:)';

Nchan = length(grad.label);
Ndip  = length(dip.inside);

if isfield(dip, 'subspace')
  % remember the original data prior to the voxel dependant subspace projection
  dat_pre_subspace = dat;
  fprintf('using subspace projection\n');
end

progress('init', feedback, 'computing inverse');
for i=1:length(dip.inside)

  progress(i/length(dip.inside), 'computing inverse %d/%d\n', i, length(dip.inside));
  i = dip.inside(i);
  
  if isfield(dip, 'leadfield')
    % reuse the leadfield that was previously computed
    lf = dip.leadfield{i};
  else
    % compute the leadfield
    lf = compute_leadfield(dip.pos(i,:), grad, vol);
  end

  if isfield(dip, 'subspace')
    % do subspace projection of the forward model
    lf = dip.subspace{i} * lf;
    % the data and the covariance become voxel dependent due to the projection
    dat = dip.subspace{i} * dat_pre_subspace;
  end
  
  % compute spatiotemporal inverse using regional source
  lfi    = pinv(lf);
  mom{i} = lfi * dat;
  rv(i)  = sum(sum((dat - lf*mom{i}).^2, 1), 2)./sum(sum(dat.^2, 1), 2);

  % for plotting convenience also compute power at each location
  % FIXME is this normalization correct?
  pow(i) = mean(sum(mom{i}(:).^2, 1));
end
progress('close');

% locations outside the head get assigned an 
for i=dip.outside
  mom{i} = [];
  rv(i)  = nan;
  pow(i) = nan;
end

% assign the output data
dipout.mom = mom(:);  % ensure that it is a column vector
dipout.rv  = rv(:);   % ensure that it is a column vector
dipout.pow = pow(:);  % ensure that it is a column vector

% add other descriptive information to the output source model
dipout.pos     = dip.pos;
dipout.inside  = dip.inside;
dipout.outside = dip.outside;

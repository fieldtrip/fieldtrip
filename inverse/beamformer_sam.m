function [dipout] = beamformer_sam(dip, sens, headmodel, dat, all_cov, varargin)

% BEAMFORMER_SAM scans on pre-defined dipole locations with a single
% dipole and returns the CTF Synthetic Aperture Magnetometry (SAM)
% beamformer estimates. location. Dipole locations that are outside
% the head will return a NaN value.
%
% Use as
%   [dipout] = beamformer_sam(dipin, sens, headmodel, dat, cov, varargin)
% where
%   dipin       is the input dipole model
%   sens        is the gradiometer definition
%   headmodel   is the volume conductor definition
%   dat         is the data matrix with the ERP or ERF
%   cov         is the data covariance or cross-spectral density matrix
% and
%   dipout      is the resulting dipole model with all details
%
% The input dipole model consists of
%   dipin.pos   positions for dipole, e.g. regular grid
%   dipin.mom   dipole orientation (optional)
%
% Additional options should be specified in key-value pairs and can be
%   ...
%
% These options influence the forward computation of the leadfield
%   'reducerank'      = reduce the leadfield rank, can be 'no' or a number (e.g. 2)
%   'normalize'       = normalize the leadfield
%   'normalizeparam'  = parameter for depth normalization (default = 0.5)

% Copyright (C) 2009,      Robert Oostenveld
% Copyright (C) 2005-2009, Arjan Hillebrand
% Copyright (C) 2005-2009, Gareth Barnes
%
% $Id$

if mod(nargin-5,2)
  % the first 5 arguments are fixed, the other arguments should come in pairs
  error('invalid number of optional arguments');
end

% get the optional input arguments
meansphereorigin  = ft_getopt(varargin, 'meansphereorigin');
feedback          = ft_getopt(varargin, 'feedback', 'text');
lambda            = ft_getopt(varargin, 'lambda', 0);
fixedori          = ft_getopt(varargin, 'fixedori', 'spinning');
% these settings pertain to the forward model, the defaults are set in ft_compute_leadfield
reducerank        = ft_getopt(varargin, 'reducerank');
normalize         = ft_getopt(varargin, 'normalize');
normalizeparam    = ft_getopt(varargin, 'normalizeparam');

% determine the mean sphere origin, required for spinning
if isempty(meansphereorigin)
  switch ft_voltype(headmodel)
    case 'singlesphere'
      meansphereorigin = headmodel.o;
    case 'localspheres'
      meansphereorigin = mean(headmodel.o, 1);
    case 'singleshell'
      meansphereorigin = mean(headmodel.bnd.pos,1);
    otherwise
      error('unsupported voltype for determining the mean sphere origin')
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dipole positions that are inside/outside the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(dip, 'inside')
  dip.inside = ft_inside_vol(dip.pos, headmodel);
end

if any(dip.inside>1)
  % convert to logical representation
  tmp = false(size(dip.pos,1),1);
  tmp(dip.inside) = true;
  dip.inside = tmp;
end

% keep the original details on inside and outside positions
originside = dip.inside;
origpos    = dip.pos;

% select only the dipole positions inside the brain for scanning
dip.pos    = dip.pos(originside,:);

if isfield(dip, 'mom')
  dip.mom = dip.mom(:, dip.inside);
end
if isfield(dip, 'leadfield')
  fprintf('using precomputed leadfields\n');
  dip.leadfield = dip.leadfield(dip.inside);
end
if isfield(dip, 'filter')
  fprintf('using precomputed filters\n');
  dip.filter = dip.filter(dip.inside);
end
dip.inside = true(size(dip.pos,1),1);

isrankdeficient = (rank(all_cov)<size(all_cov,1));

% it is difficult to give a quantitative estimate of lambda, therefore also
% support relative (percentage) measure that can be specified as string (i.e. '10%')
if ~isempty(lambda) && ischar(lambda) && lambda(end)=='%'
  ratio = sscanf(lambda, '%f%%');
  ratio = ratio/100;
  lambda = ratio * trace(all_cov)/size(all_cov,1);
end

% estimate the noise power, which is further assumed to be equal and uncorrelated over channels
if isrankdeficient
  % estimated noise floor is equal to or higher than lambda
  noise = lambda;
else
  % estimate the noise level in the covariance matrix by the smallest singular value
  noise = svd(all_cov);
  noise = noise(end);
  % estimated noise floor is equal to or higher than lambda
  noise = max(noise, lambda);
end

% the inverse only has to be computed once for all dipoles
inv_cov = pinv(all_cov + lambda * eye(size(all_cov)));
noise_cov = noise * eye(size(all_cov));

% start the scanning with the proper metric
ft_progress('init', feedback, 'scanning grid');

% the angles are the same for all dipole locations
all_angles = 0:pi/72:pi;

for diplop=1:size(dip.pos,1)
  
  vox_pos = dip.pos(diplop,:);
  
  if isfield(dip, 'leadfield')
    % reuse the leadfield that was previously computed
    lf = dip.leadfield{diplop};
  elseif isfield(dip, 'mom')
    % compute the leadfield for a fixed dipole orientation
    % FIXME this probably won't work because of the spinning/nonspinning source orientation estimate being part of this code
    lf = ft_compute_leadfield(vox_pos, sens, headmodel, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam) * dip.mom(:,diplop);
  else
    % compute the leadfield
    lf = ft_compute_leadfield(vox_pos, sens, headmodel, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam);
  end
  
  if strcmp(fixedori, 'spinning')
    % perform a non-linear search for the optimum orientation
    [tanu, tanv] = calctangent(vox_pos - meansphereorigin); % get tangential components
    % get a decent starting guess
    all_costfun_val = zeros(size(all_angles));
    for i=1:length(all_angles)
      costfun_val        = SAM_costfun(all_angles(i), vox_pos, tanu, tanv, lf, all_cov, inv_cov, noise_cov);
      all_costfun_val(i) = costfun_val;
    end
    [junk, min_ind] = min(all_costfun_val);
    
    optim_options = optimset('Display', 'final', 'TolX', 1e-3, 'Display', 'off');
    [opt_angle, fval, exitflag, output] = fminsearch(@SAM_costfun, all_angles(min_ind), optim_options, vox_pos, tanu, tanv, lf, all_cov, inv_cov, noise_cov);
    MDip        = settang(opt_angle, tanu, tanv);
    MagDip      = sqrt(dot(MDip,MDip));
    opt_vox_or  = (MDip/MagDip)';
    
    % figure
    % plot(all_angles, all_costfun_val, 'k-'); hold on; plot(opt_angle, fval, 'g*')
    % drawnow
    
  else
    % Use Sekihara's method of finding the optimum orientation
    %
    % Sekihara et al. Asymptotic SNR of scalar and vector minimum-variance
    % beamformers for neuromagnetic source reconstruction. IEEE Trans. Biomed.
    % Eng, No 10, Vol. 51, 2004, pp 1726-1734
    
    % Compute the lead field for 3 orthogonal components
    L = lf; % see above
    
    switch fixedori
      case 'gareth'
        % Compute Y1 = L' R(^-1) * L
        % and     Y2 = L' R(^-2) * L
        Y1 = L' * inv_cov * L;
        Y2 = L' * (inv_cov * inv_cov) * L;
        % find the eigenvalues and eigenvectors
        [U,S] = eig(Y2,Y1);
        
      case 'robert'
        [U,S] = svd(real(pinv(L' * inv_cov * L)));
        
      case 'stephen'
        %% Stephen Robinsons stuff? this did not work!
        L2_inv = inv(L2);
        Z2 = L2_inv * L1 * L2_inv;
        [U,S] = svds(Z2,1,0);
        % find the smallest eigenvalue and eigenvector
        Y = L2_inv*U;
        Y = Y./sqrt(dot(Y,Y));
        U = Y;
        
      otherwise
        error('unknown orimethod');
    end
    
    % The optimum orientation is the eigenvector that corresponds to the
    % smallest eigenvalue.
    
    % Double check that this is the case, because for single sphere head
    % model, one of the eigenvectors corresponds to the radial direction,
    % giving lead fields that are zero (to within machine precission).
    % The eigenvalue corresponding to this eigenvector can actually be
    % the smallest and can give the optimum (but wrong) Z-value!)
    
    ori1 = U(:,1); ori1 = ori1/norm(ori1);
    ori2 = U(:,2); ori2 = ori2/norm(ori2);
    % ori3 = U(:,3); ori3 = ori3/norm(ori3);
    
    L1 = L * ori1;
    L2 = L * ori2;
    % L3 = L * ori3;
    
    if (norm(L1)/norm(L2)) < 1e-6
      % the first orientation seems to be the silent orientation
      % use the second orientation instead
      opt_vox_or = ori2;
    else
      opt_vox_or = ori1;
    end
    
  end, % if fixedori
  
  % compute the spatial filter for the optimal source orientation
  gain        = lf * opt_vox_or;
  trgain_invC = gain' * inv_cov;
  SAMweights  = trgain_invC / (trgain_invC * gain);
  
  % remember all output details for this dipole
  dipout.pow(diplop)    = SAMweights * all_cov * SAMweights';
  dipout.noise(diplop)  = SAMweights * noise_cov * SAMweights';
  dipout.ori{diplop}    = opt_vox_or;
  dipout.filter{diplop} = SAMweights;
  if ~isempty(dat)
    dipout.mom{diplop} = SAMweights * dat;
  end
  
  ft_progress(diplop/size(dip.pos,1), 'scanning grid %d/%d\n', diplop, size(dip.pos,1));
end % for each dipole position

ft_progress('close');

% wrap it all up, prepare the complete output
dipout.inside   = originside;
dipout.pos      = origpos;

% reassign the scan values over the inside and outside grid positions
if isfield(dipout, 'leadfield')
  dipout.leadfield( originside) = dipout.leadfield;
  dipout.leadfield(~originside) = {[]};
end
if isfield(dipout, 'filter')
  dipout.filter( originside) = dipout.filter;
  dipout.filter(~originside) = {[]};
end
if isfield(dipout, 'mom')
  dipout.mom( originside) = dipout.mom;
  dipout.mom(~originside) = {[]};
end
if isfield(dipout, 'ori')
  dipout.ori( originside) = dipout.ori;
  dipout.ori(~originside) = {[]};
end
if isfield(dipout, 'pow')
  dipout.pow( originside) = dipout.pow;
  dipout.pow(~originside) = nan;
end
if isfield(dipout, 'noise')
  dipout.noise( originside) = dipout.noise;
  dipout.noise(~originside) = nan;
end

return % end of beamformer_sam() main function

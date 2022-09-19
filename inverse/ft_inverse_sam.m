function [estimate] = ft_inverse_sam(sourcemodel, sens, headmodel, dat, C, varargin)

% FT_INVERSE_SAM scans on pre-defined dipole locations with a single dipole and
% returns the Synthetic Aperture Magnetometry (SAM) beamformer estimates.
%
% Use as
%   [estimate] = ft_inverse_sam(sourcemodel, sens, headmodel, dat, cov, ...)
% where
%   sourcemodel is the input source model, see FT_PREPARE_SOURCEMODEL
%   sens        is the gradiometer or electrode definition, see FT_DATATYPE_SENS
%   headmodel   is the volume conductor definition, see FT_PREPARE_HEADMODEL
%   dat         is the data matrix with the ERP or ERF
%   cov         is the data covariance or cross-spectral density matrix
% and
%   estimate    contains the estimated source parameters
%
% Additional input arguments should be specified as key-value pairs and can include
%   'meansphereorigin'
%   'feedback'
%   'fixedori'
%   'noisecov'
%   'toi'
%
% These options influence the forward computation of the leadfield
%   'reducerank'      = 'no' or number  (default = 3 for EEG, 2 for MEG)
%   'backproject'     = 'yes' or 'no', in the case of a rank reduction this parameter determines whether the result will be backprojected onto the original subspace (default = 'yes')
%   'normalize'       = 'no', 'yes' or 'column' (default = 'no')
%   'normalizeparam'  = parameter for depth normalization (default = 0.5)
%   'weight'          = number or Nx1 vector, weight for each dipole position to compensate for the size of the corresponding patch (default = 1)
%
% See also FT_SOURCEANALYSIS, FT_PREPARE_HEADMODEL, FT_PREPARE_SOURCEMODEL

% Copyright (C) 2009,      Robert Oostenveld
% Copyright (C) 2005-2009, Arjan Hillebrand
% Copyright (C) 2005-2009, Gareth Barnes
%
% $Id$

if mod(nargin-5,2)
  % the first 5 arguments are fixed, the other arguments should come in pairs
  ft_error('invalid number of optional arguments');
end

% get the optional input arguments, or use defaults
meansphereorigin  = ft_getopt(varargin, 'meansphereorigin');
feedback          = ft_getopt(varargin, 'feedback', 'text');
fixedori          = ft_getopt(varargin, 'fixedori', 'spinning');
noisecov          = ft_getopt(varargin, 'noisecov');
toi               = ft_getopt(varargin, 'toi');

% construct the low-level options for the covariance matrix inversion as key-value pairs, these are passed to FT_INV
invopt = {};
invopt = ft_setopt(invopt, 'lambda',    ft_getopt(varargin, 'lambda', 0));
invopt = ft_setopt(invopt, 'kappa',     ft_getopt(varargin, 'kappa'));
invopt = ft_setopt(invopt, 'tolerance', ft_getopt(varargin, 'tol'));
invopt = ft_setopt(invopt, 'method',    ft_getopt(varargin, 'invmethod'));

% construct the low-level options for the leadfield computation as key-value pairs, these are passed to FT_COMPUTE_LEADFIELD
leadfieldopt = {};
leadfieldopt = ft_setopt(leadfieldopt, 'reducerank',     ft_getopt(varargin, 'reducerank'));
leadfieldopt = ft_setopt(leadfieldopt, 'backproject',    ft_getopt(varargin, 'backproject'));
leadfieldopt = ft_setopt(leadfieldopt, 'normalize',      ft_getopt(varargin, 'normalize'));
leadfieldopt = ft_setopt(leadfieldopt, 'normalizeparam', ft_getopt(varargin, 'normalizeparam'));
leadfieldopt = ft_setopt(leadfieldopt, 'weight',         ft_getopt(varargin, 'weight'));

if strcmp(fixedori, 'moiseev') && isempty(noisecov)
  ft_error('fixedori method ''moiseev'' requires a noise covariance matrix to be specified');
end

if ~strcmp(fixedori, 'spinning')
  leadfieldopt = ft_setopt(leadfieldopt, 'backproject', 'no');
end

if strcmp(fixedori, 'spinning')
  % determine the mean sphere origin, required for spinning
  % FIXME this currently does not work with precomputed leadfields, since
  % the recent updates to ft_sourceanalysis throw away the headmodel in
  % that case
  if isempty(meansphereorigin)
    switch ft_headmodeltype(headmodel)
      case 'singlesphere'
        meansphereorigin = headmodel.o;
      case 'localspheres'
        meansphereorigin = mean(headmodel.o, 1);
      case 'singleshell'
        meansphereorigin = mean(headmodel.bnd.pos,1);
      otherwise
        ft_error('unsupported headmodel type for determining the mean sphere origin')
    end
  end
end

% flags to avoid calling isfield repeatedly in the loop over grid positions (saves a lot of time)
hasmom        = isfield(sourcemodel, 'mom');
hasleadfield  = isfield(sourcemodel, 'leadfield');
hasfilter     = isfield(sourcemodel, 'filter');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dipole positions that are inside/outside the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(sourcemodel, 'inside')
  sourcemodel.inside = ft_inside_headmodel(sourcemodel.pos, headmodel);
end

if any(sourcemodel.inside>1)
  % convert to logical representation
  tmp = false(size(sourcemodel.pos,1),1);
  tmp(sourcemodel.inside) = true;
  sourcemodel.inside = tmp;
end

% keep the original details on inside and outside positions
originside = sourcemodel.inside;
origpos    = sourcemodel.pos;

% select only the dipole positions inside the brain for scanning
sourcemodel.pos = sourcemodel.pos(originside,:);

if hasmom
  ft_warning('this probably will not work because of the spinning/nonspinning source orientation estimate being part of this code');
  sourcemodel.mom = sourcemodel.mom(:, sourcemodel.inside);
end

if hasfilter
  ft_info('using precomputed filters\n');
  sourcemodel.filter = sourcemodel.filter(sourcemodel.inside);
elseif hasleadfield
  ft_info('using precomputed leadfields\n');
  sourcemodel.leadfield = sourcemodel.leadfield(sourcemodel.inside);
  
  % check that LF backprojection is not used
  lfdim  = cellfun('size', sourcemodel.leadfield, 2);
  lfrank = cellfun(@rank,  sourcemodel.leadfield);
  if ~strcmp(fixedori, 'spinning') && any(lfdim > lfrank)
    % case analytical method used, check that LF are full rank or remove it
    ft_warning('SAM with one of the analytical methods for dipole orientation estimation does not support backprojected leadfields\n');
    
    for i=1:numel(sourcemodel.leadfield)
      lf      = sourcemodel.leadfield{i};
      [U,S,V] = svd(lf, 'econ');
      sourcemodel.leadfield{i} = lf*V(:,1:2);
    end
  end
else
  ft_info('computing forward model on the fly\n');
end

if strcmp(fixedori, 'moiseev')  && ~isempty(toi)
  ft_info('computing an event-related SAM beamformer... \n')
end

isrankdeficient = (rank(C)<size(C,1));

% it is difficult to give a quantitative estimate of lambda, therefore also
% support relative (percentage) measure that can be specified as string (i.e. '10%')
lambda = ft_getopt(invopt, 'lambda');
if ~isempty(lambda) && ischar(lambda) && lambda(end)=='%'
  ratio  = sscanf(lambda, '%f%%');
  ratio  = ratio/100;
  lambda = ratio * trace(C)/size(C,1);
  invopt = ft_setopt(invopt, 'lambda', lambda);
end

% estimate the noise power, which is further assumed to be equal and uncorrelated over channels
if isrankdeficient
  % estimated noise floor is equal to or higher than lambda
  noise = lambda;
else
  % estimate the noise level in the covariance matrix by the smallest singular value
  noise = svd(C);
  noise = noise(end);
  % estimated noise floor is equal to or higher than lambda
  noise = max(noise, lambda);
end

% the inverse only has to be computed once for all dipoles
inv_cov   = ft_inv(C, invopt{:});
if isempty(noisecov)
  noisecov = noise * eye(size(C));
end


% the angles are the same for all dipole locations
all_angles = 0:pi/72:pi;

% start the scanning
ft_progress('init', feedback, 'scanning grid');
for i=1:size(sourcemodel.pos,1)
  ft_progress(i/size(sourcemodel.pos,1), 'scanning grid %d/%d\n', i, size(sourcemodel.pos,1));
  
  vox_pos = sourcemodel.pos(i,:);
  
  if hasfilter
    SAMweights = sourcemodel.filter{i};
    if size(SAMweights,1)>1
      ft_error('unsupported dimensionality of precomputed spatial filters');
    end
    
  else
    if hasleadfield && hasmom && size(sourcemodel.mom, 1)==size(sourcemodel.leadfield{i}, 2)
      % reuse the leadfield that was previously computed and project
      lf = sourcemodel.leadfield{i} * sourcemodel.mom(:,i);
    elseif  hasleadfield &&  hasmom
      % reuse the leadfield that was previously computed but don't project
      lf = sourcemodel.leadfield{i};
    elseif  hasleadfield && ~hasmom
      % reuse the leadfield that was previously computed
      lf = sourcemodel.leadfield{i};
    elseif ~hasleadfield &&  hasmom
      % compute the leadfield for a fixed dipole orientation
      lf = ft_compute_leadfield(sourcemodel.pos(i,:), sens, headmodel, leadfieldopt{:}) * sourcemodel.mom(:,i);
    else
      % compute the leadfield
      lf = ft_compute_leadfield(sourcemodel.pos(i,:), sens, headmodel, leadfieldopt{:});
    end
    
    switch fixedori
      case 'spinning'
        % perform a non-linear search for the optimum orientation
        [tanu, tanv] = calctangent(vox_pos - meansphereorigin); % get tangential components
        % get a decent starting guess
        all_costfun_val = zeros(size(all_angles));
        for j=1:length(all_angles)
          costfun_val        = SAM_costfun(all_angles(j), vox_pos, tanu, tanv, lf, C, inv_cov, noisecov);
          all_costfun_val(j) = costfun_val;
        end
        [junk, min_ind] = min(all_costfun_val);
        
        optim_options = optimset('Display', 'final', 'TolX', 1e-3, 'Display', 'off');
        [opt_angle, fval, exitflag, output] = fminsearch(@SAM_costfun, all_angles(min_ind), optim_options, vox_pos, tanu, tanv, lf, C, inv_cov, noisecov);
        MDip        = settang(opt_angle, tanu, tanv);
        MagDip      = sqrt(dot(Msourcemodel,MDip));
        opt_vox_or  = (MDip/MagDip)';
        
        % figure
        % plot(all_angles, all_costfun_val, 'k-'); hold on; plot(opt_angle, fval, 'g*')
        % drawnow
        
      case 'gareth'
        % Compute Y1 = lf' R(^-1) * lf
        % and     Y2 = lf' R(^-2) * lf
        Y1 = lf' * inv_cov * lf;
        Y2 = lf' * (inv_cov * inv_cov) * lf;
        % find the eigenvalues and eigenvectors
        [U,S] = eig(Y1,Y2);
        
      case 'robert'
        % Use Sekihara's method of finding the optimum orientation
        %
        % Sekihara et al. Asymptotic SNR of scalar and vector minimum-variance
        % beamformers for neuromagnetic source reconstruction. IEEE Trans. Biomed.
        % Eng, No 10, Vol. 51, 2004, pp 1726-1734
        [U,S] = svd(pinv(lf' * inv_cov * lf)); % JM added pinv, in order to take the orientation that belong to the maximum S
        
      case'moiseev'
        if ~isempty(toi)
          % use trial averaged variance matrix within a time of interest
          Avg   = dat(:,toi(1):toi(2));
          Cavg  = Avg*Avg'/size(Avg,2);
          Sproj = inv_cov * Cavg * inv_cov';
        else
          % in case of induced data, just use the Signal cov only
          Sproj = inv_cov;
        end
        
        Nproj = inv_cov * noisecov * inv_cov';
        
        Y1 = lf' * Sproj * lf;
        Y2 = lf' * Nproj * lf;
        % find the eigenvalues and eigenvectors
        [U,S] = eig(Y1,Y2);
        
      otherwise
        ft_error(sprintf('unknown fixedori method %s', fixedori));
    end
    
    % convert the U matrix into the optimal orientation vector
    switch fixedori
      case 'spinning'
        % do nothing, optimum orientation is already computed above
        
      otherwise
        % The optimum orientation is the eigenvector that corresponds to the
        % biggest eigenvalue (biggest value is more logical, as it relates to SNR).
        
        if isfield(headmodel,'singlesphere') && ~strcmp(fixedori,'moiseev')
          % If baseline noise covariance is not used, for single sphere head
          % model, one of the eigenvectors corresponds to the radial direction,
          % giving lead fields that are zero (to within machine precission).
          % The eigenvalue corresponding to this eigenvector can actually be
          % the biggest and can give the optimum (but wrong) Z-value!)
          ori1 = U(:,1); ori1 = ori1/norm(ori1);
          ori2 = U(:,2); ori2 = ori2/norm(ori2);
          % ori3 = U(:,3); ori3 = ori3/norm(ori3);
          
          L1 = lf * ori1;
          L2 = lf * ori2;
          
          if (norm(L1)/norm(L2)) < 1e-6
            % the first orientation seems to be the silent orientation
            % use the second orientation instead
            opt_vox_or = ori2;
          else
            opt_vox_or = ori1;
          end
          
        else
          % select eigenvector with biggest eigenvalue
          [dum, ori_inx] = sort(diag(S), 'descend');
          ori = U(:,ori_inx(1));
          opt_vox_or = ori/norm(ori);
        end
    end
    estimate.ori{i} = opt_vox_or;
    
    % compute the spatial filter for the optimal source orientation
    gain        = lf * opt_vox_or;
    trgain_invC = gain' * inv_cov;
    SAMweights  = trgain_invC / (trgain_invC * gain);
    
  end % if hasfilter or not
  
  % remember all output details for this dipole
  estimate.pow(i)    = SAMweights * C  * SAMweights';
  estimate.noise(i)  = SAMweights * noisecov * SAMweights';
  estimate.filter{i} = SAMweights;
  if ~isempty(dat)
    estimate.mom{i}  = SAMweights * dat;
  end
  if strcmp(fixedori,'moiseev') && exist('gain', 'var')
    % get pseudoZ
    Ng                  = gain' * Nproj * gain;
    Sg                  = gain' * Sproj * gain;
    estimate.pseudoZ(i) = Sg / Ng;
  end
  
end % for each dipole position
ft_progress('close');

% reassign the estimated values over the inside and outside grid positions
estimate.inside   = originside;
estimate.pos      = origpos;
if isfield(estimate, 'leadfield')
  estimate.leadfield( originside) = estimate.leadfield;
  estimate.leadfield(~originside) = {[]};
end
if isfield(estimate, 'filter')
  estimate.filter( originside) = estimate.filter;
  estimate.filter(~originside) = {[]};
end
if isfield(estimate, 'mom')
  estimate.mom( originside) = estimate.mom;
  estimate.mom(~originside) = {[]};
end
if isfield(estimate, 'ori')
  estimate.ori( originside) = estimate.ori;
  estimate.ori(~originside) = {[]};
end
if isfield(estimate, 'pow')
  estimate.pow( originside) = estimate.pow;
  estimate.pow(~originside) = nan;
end
if isfield(estimate, 'noise')
  estimate.noise( originside) = estimate.noise;
  estimate.noise(~originside) = nan;
end
if isfield(estimate, 'pseudoZ')
  estimate.pseudoZ( originside) = estimate.pseudoZ;
  estimate.pseudoZ(~originside) = nan;
end

function [dipout] = beamformer_sam(dip, sens, headmodel, dat, C, varargin)

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
  ft_error('invalid number of optional arguments');
end

% get the optional input arguments
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dipole positions that are inside/outside the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(dip, 'inside')
  dip.inside = ft_inside_headmodel(dip.pos, headmodel);
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
  ft_info('using precomputed leadfields\n');
  dip.leadfield = dip.leadfield(dip.inside);
  
  % check that LF backprojection is not used
  lfdim  = cellfun('size', dip.leadfield, 2);
  lfrank = cellfun(@rank,  dip.leadfield);
  if ~strcmp(fixedori, 'spinning') && any(lfdim > lfrank)
    % case analytical method used, check that LF are full rank or remove it
    ft_warning('SAM with one of the analytical methods for dipole orientation estimation does not support backprojected leadfields\n');
    
    for i=1:numel(dip.leadfield)
      lf      = dip.leadfield{i};
      [U,S,V] = svd(lf, 'econ');
      dip.leadfield{i} = lf*V(:,1:2);
    end
  end
end
if isfield(dip, 'filter')
  ft_info('using precomputed filters\n');
  dip.filter = dip.filter(dip.inside);
elseif strcmp(fixedori,'moiseev')  && ~isempty(toi)
  ft_info('Computing an Event Related SAM beamformer... \n')
end
dip.inside = true(size(dip.pos,1),1);

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

% start the scanning with the proper metric
ft_progress('init', feedback, 'scanning grid');

% the angles are the same for all dipole locations
all_angles = 0:pi/72:pi;

hasfilter    = isfield(dip, 'filter');
hasleadfield = isfield(dip, 'leadfield');
hasmom       = isfield(dip, 'hasmom');

for diplop=1:size(dip.pos,1)
  
  vox_pos = dip.pos(diplop,:);
  
  if hasfilter
    SAMweights = dip.filter{diplop};
    if size(SAMweights,1)>1
      ft_error('unsupported dimensionality of precomputed spatial filters');
    end
  else
    
    if hasleadfield
      % reuse the leadfield that was previously computed
      lf = dip.leadfield{diplop};
    elseif hasmom
      % compute the leadfield for a fixed dipole orientation
      % FIXME this probably won't work because of the spinning/nonspinning source orientation estimate being part of this code
      lf = ft_compute_leadfield(vox_pos, sens, headmodel, leadfieldopt{:}) * dip.mom(:,diplop);
    else
      % compute the leadfield
      lf = ft_compute_leadfield(vox_pos, sens, headmodel, leadfieldopt{:});
    end
    
    switch fixedori
      case 'spinning'
        % perform a non-linear search for the optimum orientation
        [tanu, tanv] = calctangent(vox_pos - meansphereorigin); % get tangential components
        % get a decent starting guess
        all_costfun_val = zeros(size(all_angles));
        for i=1:length(all_angles)
          costfun_val        = SAM_costfun(all_angles(i), vox_pos, tanu, tanv, lf, C, inv_cov, noisecov);
          all_costfun_val(i) = costfun_val;
        end
        [junk, min_ind] = min(all_costfun_val);
        
        optim_options = optimset('Display', 'final', 'TolX', 1e-3, 'Display', 'off');
        [opt_angle, fval, exitflag, output] = fminsearch(@SAM_costfun, all_angles(min_ind), optim_options, vox_pos, tanu, tanv, lf, C, inv_cov, noisecov);
        MDip        = settang(opt_angle, tanu, tanv);
        MagDip      = sqrt(dot(MDip,MDip));
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
          [~, ori_inx] = sort(diag(S), 'descend');
          ori = U(:,ori_inx(1));
          opt_vox_or = ori/norm(ori);
        end
    end
    dipout.ori{diplop}    = opt_vox_or;
       
    % compute the spatial filter for the optimal source orientation
    gain        = lf * opt_vox_or;
    trgain_invC = gain' * inv_cov;
    SAMweights  = trgain_invC / (trgain_invC * gain);  
    
  end
 
  % remember all output details for this dipole
  dipout.pow(diplop)    = SAMweights * C  * SAMweights';
  dipout.noise(diplop)  = SAMweights * noisecov * SAMweights'; 
  dipout.filter{diplop} = SAMweights;
  if ~isempty(dat)
    dipout.mom{diplop} = SAMweights * dat;
  end
  if strcmp(fixedori,'moiseev') && exist('gain', 'var')
    % get pseudoZ
    Ng                      = gain' * Nproj * gain;
    Sg                      = gain' * Sproj * gain; 
    dipout.pseudoZ(diplop)  = Sg / Ng;  
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
if isfield(dipout, 'pseudoZ')
  dipout.pseudoZ( originside) = dipout.pseudoZ;
  dipout.pseudoZ(~originside) = nan;
end

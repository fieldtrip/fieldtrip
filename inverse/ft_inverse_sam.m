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
%   'feedback'
%   'fixedori'          deprecated, control behaviour via 'reducerank' instead
%   'noisecov'
%   'toi'
%
% If no orientation is specified, the SAM beamformer will try to estimate the orientation from the data.
% The beamformer will either try to estimate the whole orientation, or only its tangential component.
% This is controlled by the 'reducerank' parameter. For reducerank=3, the whole orientation is estimated,
% and for reducerank=2 only the tangential component is estimated, based on an svd of the dipole's leadfield,
% treating the 3d component as the 'radial' orientation.
%
% These options influence the forward computation of the leadfield, if it has not yet been precomputed
%   'reducerank'      = 'no' or number  (default = 3 for EEG, 2 for MEG)
%   'backproject'     = 'yes' or 'no', in the case of a rank reduction this parameter determines whether the result will be backprojected onto the original subspace (default = 'yes')
%   'normalize'       = 'no', 'yes' or 'column' (default = 'no')
%   'normalizeparam'  = parameter for depth normalization (default = 0.5)
%   'weight'          = number or Nx1 vector, weight for each dipole position to compensate for the size of the corresponding patch (default = 1)
%
% See also FT_SOURCEANALYSIS, FT_PREPARE_HEADMODEL, FT_PREPARE_SOURCEMODEL

% Copyright (C) 2023,      Malte Hoeltershinken
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
feedback          = ft_getopt(varargin, 'feedback', 'text');
fixedori          = ft_getopt(varargin, 'fixedori');
noisecov          = ft_getopt(varargin, 'noisecov');
toi               = ft_getopt(varargin, 'toi');

% construct the low-level options for the covariance matrix inversion as key-value pairs, these are passed to FT_INV
invopt = {};
invopt = ft_setopt(invopt, 'lambda',    ft_getopt(varargin, 'lambda', 0));
invopt = ft_setopt(invopt, 'kappa',     ft_getopt(varargin, 'kappa'));
invopt = ft_setopt(invopt, 'tolerance', ft_getopt(varargin, 'tol'));
invopt = ft_setopt(invopt, 'method',    ft_getopt(varargin, 'invmethod'));

% backwards compatibility information
if ~isempty(fixedori)
  switch fixedori
    case 'robert'
      ft_error('fixedori=''robert'' is not supported anymore. The exact same functionality is implemented in ft_inverse_lcmv, using fixedori = ''yes''');
    case 'gareth'
      ft_warning('fixedori=''gareth'' is deprecated. Using a default noise covariance matrix, or specifying a noise covariance matrix of the form sigma * I, will reproduce the old ''gareth''-behaviour');
    case 'moiseev'
      ft_warning('fixedori=''moiseev'' is deprecated. This is now the default behaviour.');
    case 'spinning'
      ft_warning('fixedori=''spinning'' is deprecated. To estimate tangential orientations, please specify ''reducerank''=2, see the code for details.');
    otherwise
      ft_error('unknown value for fixedori option. Additionally, note that the parameter ''fixedori'' is deprecated.');
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
  ft_warning('using user-specified moments, no data-based estimation of orientation is performed');
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
  if any(lfdim > lfrank)
    % if LF was backprojected, we perform the orientation estimation only in the tangential plane
    ft_warning('at least one of the supplied leadfields does not have full rank, falling back to tangential orientation estimation\n');
    istangentialestimation = true;
  else
    istangentialestimation = false;
  end
else
  ft_info('computing forward model on the fly\n');

   % construct the low-level options for the leadfield computation as key-value pairs, these are passed to FT_COMPUTE_LEADFIELD
  leadfieldopt = {};
  leadfieldopt = ft_setopt(leadfieldopt, 'normalize',      ft_getopt(varargin, 'normalize'));
  leadfieldopt = ft_setopt(leadfieldopt, 'normalizeparam', ft_getopt(varargin, 'normalizeparam'));
  leadfieldopt = ft_setopt(leadfieldopt, 'weight',         ft_getopt(varargin, 'weight'));
  
  % check if a tangential orientation estimation is performed, i.e. if only orientations in the tangential plane are considered.
  % if this is the case, we perform the lead field rank reduction directly in  this script, since otherwise we would not
  % have access to the tangential plane defined by the full leadfield
  reducerank = ft_getopt(varargin, 'reducerank');
  if isempty(reducerank)
    if ft_senstype(sens, 'meg')
      istangentialestimation = true;
    elseif ft_senstype(sens, 'eeg')
      istangentialestimation = false;
    else
      ft_error('unknown sensor type');
    end
  else
    if (ischar(reducerank) && strcmp(reducerank, 'yes')) || (isnumeric(reducerank) && reducerank == 2)
      istangentialestimation = true;
    elseif (ischar(reducerank) && strcmp(reducerank, 'no')) || (isnumeric(reducerank) && reducerank == 3)
      istangentialestimation = false;
    else
      ft_error('invalid option for reducerank, please specify one from {''yes'', ''no'', 2, 3}')
    end
  end % if isempty(reducerank) or not
  leadfieldopt = ft_setopt(leadfieldopt, 'reducerank', 'no');
  
  % for singlesphere MEG headmodels, the leadfield of a radial dipole is zero. Hence it does not make sense to try to estimate radial orientation components in this case.
  if strcmp(headmodel.type,'singlesphere') && ft_senstype(sens, 'meg') && ~istangentialestimation
    ft_warning('singlesphere MEG approach without rank reduction specified. To avoid numerical difficulties, we override this.');
    istangentialestimation = true;
  end
end % if hasfilter or not

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
invC = ft_inv(C, invopt{:});

% If no orientation is specified, it is estimated from the data. This function implements two approaches
% for this estimation, both based on choosing an orientation that maximizes some form of "pseudo-Z". Concretely, these two approaches maximize
% 1): pseudo-Z = projected_signal_power / projected_noise_power = (w' C w) / (w' N w),
% 2): event-related pseudo-Z = projected_evoked_activity_power / projected_noise_power = (w' Cavg w) / (w' N w), 
% where w is the spatial filter, C is the covariance matrix, N is the noise covarianc matrix, and Cavg is the
% second moment matrix of the averaged data over some time interval. For a detailled discussion
% regarding these approaches, we refer to 
% Moiseev et al., Application of multi-source minimum variance beamformers for reconstruction of correlated neural activity,
% NeuroImage, Volume 58, Issue 2, 15 September 2011, Pages 481-496.
% Concretely, 1) is a 1-dimensional version of the "MPZ" approach in that paper, and 2) is a 1-dimensional version of the "MER" approach
% Note that w is of the form c * invC * l, where l is the leadfield vector and c is some scalar constant (and hence cancels in the quotients above).
% Since the matrices C, Cavg, invC, and N do not depend on the source position, we can precompute the corresponding expressions here. 
if ~isempty(toi)
  % case 2) above
  ft_info('computing an event-related SAM beamformer... \n')
  Avg = dat(:,toi(1):toi(2));
  Cavg  = Avg*Avg'/size(Avg,2);
  signal_congruence = invC' * Cavg * invC;
else
  % case 1) above
  signal_congruence = invC;
end

if isempty(noisecov)
  ft_info('no noise covariance matrix supplied, using an estimate of the form noise_level * I'); 
  noisecov = noise * eye(size(C));
end
noise_congruence = invC' * noisecov * invC;

% start the scanning
ft_progress('init', feedback, 'scanning grid');
for i=1:size(sourcemodel.pos,1)
  ft_progress(i/size(sourcemodel.pos,1), 'scanning grid %d/%d\n', i, size(sourcemodel.pos,1));
  
  if hasfilter
    filt = sourcemodel.filter{i};
    if size(filt,1)>1
      ft_error('unsupported dimensionality of precomputed spatial filters');
    end
    
  else
    %  we want to compute the spatial filter. For this, we need a leadfield vector. There are 3 options for this.
    % 1): The leadfield vector is already specified.
    % 2): A dipole moment is specified. The leadfield vector can then be optained via multiplication with the full n x 3 - leadfield.
    % 3): Neither a leadfield vector nor a moment are specified. 
    % In the case of 3), we derive the orientation from the data by choosing the orientation that maximizes the pseudo-Z score.
    if hasleadfield && size(sourcemodel.leadfield{i}, 2) == 1
      gain = sourcemodel.leadfield{i};
    elseif hasmom
      % try to project moment to leadfield vector
      if hasleadfield
        lf = sourcemodel.leadfield{i};
      else
        lf = ft_compute_leadfield(sourcemodel.pos(i,:), sens, headmodel, leadfieldopt{:});
      end
      
      if size(sourcemodel.mom, 1) == size(lf, 2)
        gain = lf * sourcemodel.mom(:, i);
      else
        ft_error('first dimension of field "mom" does not match second dimension of leadfield');
      end
    else
      % estimate orientation from data
      if hasleadfield
        lf = sourcemodel.leadfield{i};
      else
        lf = ft_compute_leadfield(sourcemodel.pos(i,:), sens, headmodel, leadfieldopt{:});
      end

      if istangentialestimation
        [dum_U, dum_S, V_lf] = svd(lf, 'econ');
        lf = lf * V_lf(:, 1:2);
      end

      projected_signal = lf' * signal_congruence * lf;
      projected_noise = lf' * noise_congruence * lf;
      [U, S] = eig(projected_signal, projected_noise);
      [dum, ori_inx] = sort(diag(S), 'descend');

      max_ori = U(:, ori_inx(1));
      max_ori = max_ori / norm(max_ori);

      gain = lf * max_ori;

      if istangentialestimation
        estimate.ori{i} = V_lf(:, 1:2) * max_ori;
      else
        estimate.ori{i} = max_ori;
      end
    end % if hasleadfield or not

    % compute the spatial filter
    trgain_invC = gain' * invC;
    filt  = trgain_invC / (trgain_invC * gain);
    
  end % if hasfilter or not
  
  % remember all output details for this dipole
  estimate.pow(i)    = filt * C  * filt';
  estimate.noise(i)  = filt * noisecov * filt';
  estimate.filter{i} = filt;
  if ~isempty(dat)
    estimate.mom{i}  = filt * dat;
  end
  if ~isempty(toi)
    estimate.pseudoZ(i) = (filt * Cavg * filt') / estimate.noise(i);
  else
    estimate.pseudoZ(i) = estimate.pow(i) / estimate.noise(i);
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

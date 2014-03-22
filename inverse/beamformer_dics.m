function [dipout] = beamformer_dics(dip, grad, vol, dat, Cf, varargin)

% BEAMFORMER_DICS scans on pre-defined dipole locations with a single dipole
% and returns the beamformer spatial filter output for a dipole on every
% location.  Dipole locations that are outside the head will return a
% NaN value.
%
% Use as
%   [dipout] = beamformer_dics(dipin, grad, vol, dat, cov, varargin)
% where
%   dipin       is the input dipole model
%   grad        is the gradiometer definition
%   vol         is the volume conductor definition
%   dat         is the data matrix with the ERP or ERF
%   cov         is the data covariance or cross-spectral density matrix
% and
%   dipout      is the resulting dipole model with all details
%
% The input dipole model consists of
%   dipin.pos   positions for dipole, e.g. regular grid, Npositions x 3
%   dipin.mom   dipole orientation (optional), 3 x Npositions
% and can additionally contain things like a precomputed filter.
%
% Additional options should be specified in key-value pairs and can be
%  'Pr'               = power of the external reference channel
%  'Cr'               = cross spectral density between all data channels and the external reference channel
%  'refdip'           = location of dipole(s) with which coherence is computed
%                       computed. use 'all' to to treat all source dipoles as ref dipoles
%  'lambda'           = regularisation parameter
%  'powmethod'        = can be 'trace' or 'lambda1'
%  'feedback'         = give ft_progress indication, can be 'text', 'gui' or 'none'
%  'fixedori'         = use fixed or free orientation,                 can be 'yes' or 'no'
%  'projectnoise'     = project noise estimate through filter,         can be 'yes' or 'no'
%  'realfilter'       = construct a real-valued filter,                can be 'yes' or 'no'
%  'keepfilter'       = remember the beamformer filter,                can be 'yes' or 'no'
%  'keepleadfield'    = remember the forward computation,              can be 'yes' or 'no'
%  'keepcsd'          = remember the estimated cross-spectral density, can be 'yes' or 'no'
%  'connsum'          = Summary of coherence values (if more than one refdips specifled. 
%                       Can be 'max', 'absmax','sse','var','std'.
%  'imagcoh'          = Compute d
%
% These options influence the forward computation of the leadfield
%  'reducerank'       = reduce the leadfield rank, can be 'no' or a number (e.g. 2)
%  'normalize'        = normalize the leadfield
%  'normalizeparam'   = parameter for depth normalization (default = 0.5)
%
% If the dipole definition only specifies the dipole location, a rotating
% dipole (regional source) is assumed on each location. If a dipole moment
% is specified, its orientation will be used and only the strength will
% be fitted to the data.

% Copyright (C) 2003-2008, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

if mod(nargin-5,2)
  % the first 5 arguments are fixed, the other arguments should come in pairs
  error('invalid number of optional arguments');
end

% these optional settings do not have defaults
Pr             = keyval('Pr',            varargin);
Cr             = keyval('Cr',            varargin);
refdip         = keyval('refdip',        varargin);
powmethod      = keyval('powmethod',     varargin); % the default for this is set below
realfilter     = keyval('realfilter',    varargin); % the default for this is set below
% these settings pertain to the forward model, the defaults are set in compute_leadfield
reducerank     = keyval('reducerank',     varargin);
normalize      = keyval('normalize',      varargin);
normalizeparam = keyval('normalizeparam', varargin);
% these optional settings have defaults
feedback       = keyval('feedback',      varargin); if isempty(feedback),      feedback = 'text';            end
keepcsd        = keyval('keepcsd',       varargin); if isempty(keepcsd),       keepcsd = 'no';               end
keepfilter     = keyval('keepfilter',    varargin); if isempty(keepfilter),    keepfilter = 'no';            end
keepleadfield  = keyval('keepleadfield', varargin); if isempty(keepleadfield), keepleadfield = 'no';         end
lambda         = keyval('lambda',        varargin); if isempty(lambda  ),      lambda = 0;                   end
projectnoise   = keyval('projectnoise',  varargin); if isempty(projectnoise),  projectnoise = 'yes';         end
fixedori       = keyval('fixedori',      varargin); if isempty(fixedori),      fixedori = 'no';              end
subspace       = keyval('subspace',      varargin); 
connsum        = keyval('connsum',       varargin); if isempty(connsum),      connsum = 'absmax';              end
imagcoh        = keyval('imagcoh',       varargin); if isempty(imagcoh),      imagcoh = 0;              end
% convert the yes/no arguments to the corresponding logical values
keepcsd        = strcmp(keepcsd,       'yes');
keepfilter     = strcmp(keepfilter,    'yes');
keepleadfield  = strcmp(keepleadfield, 'yes');
projectnoise   = strcmp(projectnoise,  'yes');
fixedori       = strcmp(fixedori,      'yes');
dofeedback     = ~strcmp(feedback,     'none');
% FIXME besides regular/complex lambda1, also implement a real version

% default is to use the largest singular value of the csd matrix, see Gross 2001
if isempty(powmethod)
  powmethod = 'lambda1';
end

% default is to be consistent with the original description of DICS in Gross 2001
if isempty(realfilter)
  realfilter = 'no';
end

% use these two logical flags instead of doing the string comparisons each time again
powtrace   = strcmp(powmethod, 'trace');
powlambda1 = strcmp(powmethod, 'lambda1');

if ~isempty(Cr)
  % ensure that the cross-spectral density with the reference signal is a column matrix
  Cr = Cr(:);
end

if isfield(dip, 'mom') && fixedori
  error('you cannot specify a dipole orientation and fixedmom simultaneously');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dipole positions that are inside/outside the brain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(dip, 'inside') && ~isfield(dip, 'outside');
  insideLogical = ft_inside_vol(dip.pos, vol);
  dip.inside = find(insideLogical);
  dip.outside = find(~dip.inside);
elseif isfield(dip, 'inside') && ~isfield(dip, 'outside');
  dip.outside    = setdiff(1:size(dip.pos,1), dip.inside);
elseif ~isfield(dip, 'inside') && isfield(dip, 'outside');
  dip.inside     = setdiff(1:size(dip.pos,1), dip.outside);
end

% flags to avoid calling isfield repeatedly in the loop over grid positions
% (saves a lot of time)
hasmom = 0;
hasleadfield = 0;
hasfilter = 0;
hassubspace = 0;

% select only the dipole positions inside the brain for scanning
dip.origpos     = dip.pos;
dip.originside  = dip.inside;
dip.origoutside = dip.outside;
if hasmom
  hasmom = 1;
  dip.mom = dip.mom(:,dip.inside);
end
if isfield(dip, 'leadfield')
  hasleadfield = 1;
  if dofeedback
    fprintf('using precomputed leadfields\n');
  end
  dip.leadfield = dip.leadfield(dip.inside);
end
if isfield(dip, 'filter')
  hasfilter = 1;
  if dofeedback
    fprintf('using precomputed filters\n');
  end
  dip.filter = dip.filter(dip.inside);
end
if isfield(dip, 'subspace')
  hassubspace = 1;
  if dofeedback
    fprintf('using subspace projection\n');
  end
  dip.subspace = dip.subspace(dip.inside);
end
dip.pos     = dip.pos(dip.inside, :);
dip.inside  = 1:size(dip.pos,1);
dip.outside = [];

% dics has the following sub-methods, which depend on the function input arguments
% power only, cortico-muscular coherence and cortico-cortical coherence
if ~isempty(Cr) && ~isempty(Pr) && isempty(refdip)
  % compute cortico-muscular coherence, using reference cross spectral density
  submethod = 'dics_refchan';
elseif isempty(Cr) && isempty(Pr) && ~isempty(refdip)
  % compute cortico-cortical coherence with a dipole at the reference position
  submethod = 'dics_refdip';
elseif isempty(Cr) && isempty(Pr) && isempty(refdip)
  % only compute power of a dipole at the grid positions
  submethod = 'dics_power';
else
  error('invalid combination of input arguments for dics');
end

isrankdeficient = (rank(Cf)<size(Cf,1));

% it is difficult to give a quantitative estimate of lambda, therefore also
% support relative (percentage) measure that can be specified as string (e.g. '10%')
if ~isempty(lambda) && ischar(lambda) && lambda(end)=='%'
  ratio = sscanf(lambda, '%f%%');
  ratio = ratio/100;
  lambda = ratio * trace(Cf)/size(Cf,1);
end

if projectnoise
  % estimate the noise power, which is further assumed to be equal and uncorrelated over channels
  if isrankdeficient
    % estimated noise floor is equal to or higher than lambda
    noise = lambda;
  else
    % estimate the noise level in the covariance matrix by the smallest singular value
    noise = svd(Cf);
    noise = noise(end);
    % estimated noise floor is equal to or higher than lambda
    noise = max(noise, lambda);
  end
end

% the inverse only has to be computed once for all dipoles
if strcmp(realfilter, 'yes')
  % the filter is computed using only the leadfield and the inverse covariance or CSD matrix
  % therefore using the real-valued part of the CSD matrix here ensures a real-valued filter
  invCf = pinv(real(Cf) + lambda * eye(size(Cf)));
else
  invCf = pinv(Cf + lambda * eye(size(Cf)));
end

if hassubspace
  if dofeedback
    fprintf('using source-specific subspace projection\n');
  end
  % remember the original data prior to the voxel dependent subspace projection
  dat_pre_subspace = dat;
  Cf_pre_subspace = Cf;
  if strcmp(submethod, 'dics_refchan')
    Cr_pre_subspace = Cr;
    Pr_pre_subspace = Pr;
  end
elseif ~isempty(subspace)
  if dofeedback
    fprintf('using data-specific subspace projection\n');
  end
  % TODO implement an "eigenspace beamformer" as described in Sekihara et al. 2002 in HBM
  if numel(subspace)==1,
    % interpret this as a truncation of the eigenvalue-spectrum
    % if <1 it is a fraction of the largest eigenvalue
    % if >=1 it is the number of largest eigenvalues
    dat_pre_subspace = dat;
    Cf_pre_subspace  = Cf;
    [u, s, v] = svd(real(Cf));
    if subspace<1,
      sel      = find(diag(s)./s(1,1) > subspace);
      subspace = max(sel);
    end
    
    Cf       = s(1:subspace,1:subspace);
    % this is equivalent to subspace*Cf*subspace' but behaves well numerically
    % by construction.
    invCf    = diag(1./diag(Cf));
    subspace = u(:,1:subspace)';
    dat      = subspace*dat;

    if strcmp(submethod, 'dics_refchan')
      Cr = subspace*Cr;
    end
    
  else
    Cf_pre_subspace  = Cf;
    Cf    = subspace*Cf*subspace'; % here the subspace can be different from
    % the singular vectors of Cy, so we have to do the sandwiching as opposed
    % to line 216
    if strcmp(realfilter, 'yes')
      invCf = pinv(real(Cf));
    else
      invCf = pinv(Cf);
    end
  
    if strcmp(submethod, 'dics_refchan')
      Cr = subspace*Cr;
    end
  end
end

% start the scanning with the proper metric
if dofeedback
  ft_progress('init', feedback, 'scanning grid');
end
switch submethod

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dics_power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'dics_power'
    % only compute power of a dipole at the grid positions
    for i=1:size(dip.pos,1)
      if hasleadfield && hasmom && size(dip.mom, 1)==size(dip.leadfield{i}, 2)
        % reuse the leadfield that was previously computed and project
        lf = dip.leadfield{i} * dip.mom(:,i);
      elseif  hasleadfield &&  hasmom
        % reuse the leadfield that was previously computed but don't project
        lf = dip.leadfield{i};
      elseif  hasleadfield && ~hasmom
        % reuse the leadfield that was previously computed
        lf = dip.leadfield{i};    
      elseif ~hasleadfield && hasmom
        % compute the leadfield for a fixed dipole orientation
        lf = ft_compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam) * dip.mom(:,i);
      else
        % compute the leadfield
        lf = ft_compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize, 'normalizeparam', normalizeparam);
      end
      if hassubspace
        % do subspace projection of the forward model
        lf = dip.subspace{i} * lf;
        % the cross-spectral density becomes voxel dependent due to the projection
        Cf    = dip.subspace{i} * Cf_pre_subspace * dip.subspace{i}';
        if strcmp(realfilter, 'yes')
          invCf = pinv(dip.subspace{i} * (real(Cf_pre_subspace) + lambda * eye(size(Cf_pre_subspace))) * dip.subspace{i}');
        else
          invCf = pinv(dip.subspace{i} * (Cf_pre_subspace + lambda * eye(size(Cf_pre_subspace))) * dip.subspace{i}');
        end
      elseif ~isempty(subspace)
        % do subspace projection of the forward model only
        lforig = lf;
        lf     = subspace * lf;
    
        % according to Kensuke's paper, the eigenspace bf boils down to projecting
        % the 'traditional' filter onto the subspace
        % spanned by the first k eigenvectors [u,s,v] = svd(Cy); filt = ESES*filt; 
        % ESES = u(:,1:k)*u(:,1:k)';
        % however, even though it seems that the shape of the filter is identical to
        % the shape it is obtained with the following code, the w*lf=I does not
        % hold.
      end
      
      if hasfilter
        % use precomputed filter
        filt = dip.filter{i};
      else
        % compute filter
        filt = pinv(lf' * invCf * lf) * lf' * invCf;              % Gross eqn. 3, use PINV/SVD to cover rank deficient leadfield
      end
      
      if fixedori
        % use single dipole orientation
        if hasfilter && size(filt,1) == 1
          % provided precomputed filter already projects to one
          % orientation, nothing to be done here
        else
          % find out the optimal dipole orientation
          [u, s, v] = svd(real(filt * Cf * ctranspose(filt)));
          maxpowori = u(:,1);
          eta = s(1,1)./s(2,2);
          
          % and compute the leadfield for that orientation
          lf  = lf * maxpowori;
          dipout.ori{i} = maxpowori;
          dipout.eta{i} = eta;
          if ~isempty(subspace), lforig = lforig * maxpowori; end
          
          % recompute the filter to only use that orientation
          filt = pinv(lf' * invCf * lf) * lf' * invCf;
        end
      elseif hasfilter && size(filt,1) == 1
        error('the precomputed filter you provided projects to a single dipole orientation, but you request fixedori=''no''; this is invalid. Either provide a filter with the three orientations retained, or specify fixedori=''yes''.');
      end
      
      csd = filt * Cf * ctranspose(filt);                         % Gross eqn. 4 and 5
      if powlambda1
        if size(csd,1) == 1
          % only 1 orientation, no need to do svd
          dipout.pow(i) = real(csd);
        else
          dipout.pow(i) = lambda1(csd);                           % compute the power at the dipole location, Gross eqn. 8
        end
      elseif powtrace
        dipout.pow(i) = real(trace(csd));                         % compute the power at the dipole location
      end
      if keepcsd
        dipout.csd{i} = csd;
      end
      if projectnoise
        if powlambda1
          dipout.noise(i) = noise * lambda1(filt * ctranspose(filt));
        elseif powtrace
          dipout.noise(i) = noise * real(trace(filt * ctranspose(filt)));
        end
        if keepcsd
          dipout.noisecsd{i} = noise * filt * ctranspose(filt);
        end
      end
      if keepfilter
        if ~isempty(subspace)
          dipout.filter{i} = filt*subspace; %FIXME should this be subspace, or pinv(subspace)?
        else
          dipout.filter{i} = filt;
        end
      end
      if keepleadfield
        if ~isempty(subspace)
          dipout.leadfield{i} = lforig;
        else
          dipout.leadfield{i} = lf;
        end
      end
      if dofeedback
        ft_progress(i/size(dip.pos,1), 'scanning grid %d/%d\n', i, size(dip.pos,1));
      end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dics_refchan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'dics_refchan'
    % compute cortico-muscular coherence, using reference cross spectral density
    for i=1:size(dip.pos,1)
      if hasleadfield
        % reuse the leadfield that was previously computed
        lf = dip.leadfield{i};
      elseif hasmom
        % compute the leadfield for a fixed dipole orientation
        lf = ft_compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize) .* dip.mom(i,:)';
      else
        % compute the leadfield
        lf = ft_compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize);
      end
      if hassubspace
        % do subspace projection of the forward model
        lforig = lf;
        lf = dip.subspace{i} * lf;
        % the cross-spectral density becomes voxel dependent due to the projection
        Cf    = dip.subspace{i} * Cf_pre_subspace * dip.subspace{i}';
        invCf = pinv(dip.subspace{i} * (Cf_pre_subspace + lambda * eye(size(Cf))) * dip.subspace{i}');
      elseif ~isempty(subspace)
        % do subspace projection of the forward model only
        lforig = lf;
        lf     = subspace * lf;
    
        % according to Kensuke's paper, the eigenspace bf boils down to projecting
        % the 'traditional' filter onto the subspace
        % spanned by the first k eigenvectors [u,s,v] = svd(Cy); filt = ESES*filt; 
        % ESES = u(:,1:k)*u(:,1:k)';
        % however, even though it seems that the shape of the filter is identical to
        % the shape it is obtained with the following code, the w*lf=I does not
        % hold.
      end
      
      if hasfilter
        % use precomputed filter
        filt = dip.filter{i};
      else
        % compute filter
        filt = pinv(lf' * invCf * lf) * lf' * invCf;              % Gross eqn. 3, use PINV/SVD to cover rank deficient leadfield
      end
      
      if fixedori
        % use single dipole orientation
        if hasfilter && size(filt,1) == 1
          % provided precomputed filter already projects to one
          % orientation, nothing to be done here
        else
          % find out the optimal dipole orientation
          [u, s, v] = svd(real(filt * Cf * ctranspose(filt)));
          maxpowori = u(:,1);
          
          % compute the leadfield for that orientation
          lf  = lf * maxpowori;
          dipout.ori{i} = maxpowori;
          
          % recompute the filter to only use that orientation
          filt = pinv(lf' * invCf * lf) * lf' * invCf;
        end
      elseif hasfilter && size(filt,1) == 1
        error('the precomputed filter you provided projects to a single dipole orientation, but you request fixedori=''no''; this is invalid. Either provide a filter with the three orientations retained, or specify fixedori=''yes''.');
      end

      if powlambda1
        [pow, ori] = lambda1(filt * Cf * ctranspose(filt));            % compute the power and orientation at the dipole location, Gross eqn. 4, 5 and 8
      elseif powtrace
        pow = real(trace(filt * Cf * ctranspose(filt)));               % compute the power at the dipole location
      end
      csd = filt*Cr;                                                   % Gross eqn. 6
      if powlambda1
        % FIXME this should use the dipole orientation with maximum power
        coh = lambda1(csd)^2 / (pow * Pr);                             % Gross eqn. 9
      elseif powtrace
        coh = norm(csd)^2 / (pow * Pr);
      end
      dipout.pow(i) = pow;
      dipout.coh(i) = coh;
      if keepcsd
        dipout.csd{i} = csd;
      end
      if projectnoise
        if powlambda1
          dipout.noise(i) = noise * lambda1(filt * ctranspose(filt));
        elseif powtrace
          dipout.noise(i) = noise * real(trace(filt * ctranspose(filt)));
        end
        if keepcsd
          dipout.noisecsd{i} = noise * filt * ctranspose(filt);
        end
      end
      if keepfilter
        dipout.filter{i} = filt;
      end
      if keepleadfield
        if ~isempty(subspace)
          dipout.leadfield{i} = lforig;
        else
          dipout.leadfield{i} = lf;
        end
      end
      if dofeedback
        ft_progress(i/size(dip.pos,1), 'scanning grid %d/%d\n', i, size(dip.pos,1));
      end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dics_refdip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'dics_refdip'
    if hassubspace || ~isempty(subspace)
      error('subspace projections are not supported for beaming cortico-cortical coherence');
    end
    if fixedori
      error('fixed orientations are not supported for beaming cortico-cortical coherence');
    end
    % compute cortio-cortical coherence with a dipole at the reference position


% select all source positions as refereces (this does whole-brain
% connectivity)
    if isequal(refdip,'all')
        refdip=dip.pos;
    end
    
    % identify refdips that match existing source positions to reuse later
        [lf_keep keep_idx]=ismember(refdip,dip.pos,'rows');
  
        % check if any refdipoles are not part of original source spaces and compute New leadfileds for new source positions. 
        new_idx=zeros(size(keep_idx));
        lf_new={}; 
            for i=1:length(lf_keep)
                if lf_keep(i)
                 %   filt{keep_idx(i)} = pinv(lf{i}' * invCf * lf{i}) * lf{i}' * invCf;
                else
                    lf_new{end+1} = ft_compute_leadfield(refdip(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize);
               %     filt_new{end+1} = pinv(lf_new{end}' * invCf * lf_new{end}) * lf_new{end}' * invCf;
                    new_idx(i) = length(lf_new);
                end
            end
            
            % get indices for source positions thst match refdips
            [lf_keep2 keep_idx2]=ismember(dip.pos,refdip,'rows');

            
            


    dipout.coh = NaN.*ones(size(dip.pos,1),size(refdip,1));
    iscsd = ones(size(dip.pos,1),size(refdip,1));
    
               if keepcsd
                      dipout.csd = cell(size(dip.pos,1),size(refdip,1));
                  end
                  if projectnoise && keepcsd
                      dipout.noisecsd = cell(size(dip.pos,1),size(refdip,1));
                  end
                  
                  
    
    
    % construct filters for source dipoles
    for i=1:size(dip.pos,1)
      if hasleadfield
        % reuse the leadfield that was previously computed
        lf{i} = dip.leadfield{i};
      elseif hasmom
        % compute the leadfield for a fixed dipole orientation
        lf{i} = ft_compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize) .* dip.mom(i,:)';
      else
        % compute the leadfield
        lf{i} = ft_compute_leadfield(dip.pos(i,:), grad, vol, 'reducerank', reducerank, 'normalize', normalize);
      end
      if hasfilter
        % use the provided filter
        filt{i} = dip.filter{i};
      else
        % construct the spatial filter for the second dipole location
        filt{i} = pinv(lf{i}' * invCf * lf{i}) * lf{i}' * invCf;   %  use PINV/SVD to cover rank deficient leadfield
      end
      

      

    
    end
    
      
    % reconstruct csds between source and fer dipoles
    for i=1:size(dip.pos,1)  
        
              %reconstruct power for sources
      if powlambda1
        pow = lambda1(filt{i} * Cf * ctranspose(filt{i}));     % compute the power at the second dipole location, Gross eqn. 8
      elseif powtrace
        pow = real(trace(filt{i} * Cf * ctranspose(filt{i}))); % compute the power at the second dipole location
      end
      dipout.pow(i) = pow;
      
      if projectnoise
          if powlambda1
              dipout.noise(i) = noise * lambda1(filt{i} * ctranspose(filt{i}));
          elseif powtrace
              dipout.noise(i) = noise * real(trace(filt{i} * ctranspose(filt{i})));
          end
      end
      

      
      % construct filters for referece dipoles if needed
 
      for j=1:size(refdip,1)
          

          
          if new_idx(j)>0
              lf2 = lf_new{new_idx(j)};
              filt2 = lf2' * invCf * lf2 * lf2' * invCf;   %  use PINV/SVD to cover rank deficient leadfield


              
          else
              
                    if keep_idx(j)==i;
               iscsd(i,j)=NaN;
                    end
      
             % if previously computed for alternate pair , just take the
          % transpose of that
          
              if keep_idx2(i)>0
              if ~isnan(dipout.coh(keep_idx(j),keep_idx2(i)))
                  dipout.coh(i,j)=(dipout.coh(keep_idx(j),keep_idx2(i)));
                  if keepcsd
                      dipout.csd{i,j} = ctranspose(dipout.csd{keep_idx(j),keep_idx2(i)});
                  end
                  if projectnoise && keepcsd
                      dipout.noisecsd{i,j} = ctranspose(dipout.noisecsd{keep_idx(j),keep_idx2(i)});
                  end
                  if imagcoh
                  dipout.icoh(i,j)=ctranspose(dipout.icoh(keep_idx(j),keep_idx2(i)));
                  end
                  continue
              end
              
              end

              filt2 = filt{keep_idx(j)};
              lf2 = dip.leadfield{keep_idx(j)};
              

          end

          
          if powlambda1
              Pref = lambda1(filt2 * Cf * ctranspose(filt2));      % compute the power at the first dipole location, Gross eqn. 8
          elseif powtrace
              Pref = real(trace(filt2 * Cf * ctranspose(filt2)));  % compute the power at the first dipole location
          end
          

    
      
          
      csd = filt{i} * Cf * ctranspose(filt2);                % compute the cross spectral density between the two dipoles, Gross eqn. 4
      
%       if ~isempty(complexcoh) % compute imag, real or abs of csd before computing coherence
%           csd=eval([complexcoh '(csd);']);
%       end
      
      if powlambda1
        coh = lambda1(csd)^2 / (pow * Pref);               % compute the coherence between the first and second dipole
     if imagcoh
         icoh = lambda1(imag(csd))^2 / (pow * Pref);
     end
      elseif powtrace
        coh = real(trace((csd)))^2 / (pow * Pref);         % compute the coherence between the first and second dipole
     if imagcoh
        icoh = real(trace(imag(csd)))^2 / (pow * Pref);         % compute the coherence between the first and second dipole
     end
      end

      dipout.coh(i,j) = coh;
      
            if imagcoh
          dipout.icoh(i,j) = icoh;
      end
      
      if keepcsd
          dipout.csd{i,j} = csd;
      end
      if projectnoise && keepcsd
          dipout.noisecsd{i,j} = noise * filt{i} * ctranspose(filt2);
      end
      end
      


    if keepleadfield
        dipout.leadfield{i} = lf{:};
    end
    
    if keepfilter
        dipout.filter{i} = filt{:};
    end


    
    if dofeedback
        ft_progress(i/size(dip.pos,1), 'scanning grid %d/%d\n', i, size(dip.pos,1));
    end
    end
    
    switch connsum
        case 'max'
            dipout.cohsum=nanmax(dipout.coh.*iscsd,[],2);
        case 'absmax'
            dipout.cohsum=nanmax(abs(dipout.coh).*iscsd,[],2);
        case 'sumsqr'
            dipout.cohsum=nansum((dipout.coh.^2).*iscsd,[],2);
        case 'std'
            dipout.cohsum=nanstd((dipout.coh).*iscsd,[],2);
        case 'var'
            dipout.cohsum=nanstd((cdipout.cohoh).*iscsd,[],2).^2;
    end
    
    if imagcoh
            switch connsum
        case 'max'
            dipout.icohsum=nanmax(dipout.icoh.*iscsd,[],2);
        case 'absmax'
            dipout.icohsum=nanmax(abs(dipout.icoh).*iscsd,[],2);
        case 'sumsqr'
            dipout.icohsum=nansum((dipout.icoh.^2).*iscsd,[],2);
        case 'std'
            dipout.icohsum=nanstd((dipout.icoh).*iscsd,[],2);
        case 'var'
            dipout.icohsum=nanstd((dipout.icoh).*iscsd,[],2).^2;
            end
    end

    
end % switch submethod

if dofeedback
  ft_progress('close');
end

dipout.inside  = dip.originside;
dipout.outside = dip.origoutside;
dipout.pos     = dip.origpos;

% reassign the scan values over the inside and outside grid positions
if isfield(dipout, 'leadfield')
  dipout.leadfield(dipout.inside)  = dipout.leadfield;
  dipout.leadfield(dipout.outside) = {[]};
end
if isfield(dipout, 'filter')
  dipout.filter(dipout.inside)  = dipout.filter;
  dipout.filter(dipout.outside) = {[]};
end
if isfield(dipout, 'ori')
  dipout.ori(dipout.inside)  = dipout.ori;
  dipout.ori(dipout.outside) = {[]};
end
if isfield(dipout, 'pow')
  dipout.pow(dipout.inside)  = dipout.pow;
  dipout.pow(dipout.outside) = nan;
end
if isfield(dipout, 'noise')
  dipout.noise(dipout.inside)  = dipout.noise;
  dipout.noise(dipout.outside) = nan;
end
if isfield(dipout, 'coh')
  dipout.coh(dipout.inside,:)  = dipout.coh;
  dipout.coh(dipout.outside,:) = nan;
end
if isfield(dipout, 'csd')
  dipout.csd(dipout.inside,:)  = dipout.csd;
  dipout.csd(dipout.outside,:) = {[]};
end
if isfield(dipout, 'icoh')
  dipout.icoh(dipout.inside,:)  = dipout.icoh;
  dipout.icoh(dipout.outside,:) = nan;
end
if isfield(dipout, 'csdnoise')
  dipout.csd(dipout.inside,:)  = dipout.csdnoise;
  dipout.csd(dipout.outside,:) = {[]};
end
if isfield(dipout, 'cohsum')
  dipout.cohsum(dipout.inside,:)  = dipout.cohsum;
  dipout.cohsum(dipout.outside,:) = NaN;
end
if isfield(dipout, 'icohsum')
  dipout.icohsum(dipout.inside,:)  = dipout.icohsum;
  dipout.icohsum(dipout.outside,:) = NaN;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to obtain the largest singular value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s, ori] = lambda1(x)
% determine the largest singular value, which corresponds to the power along the dominant direction
[u, s, v] = svd(x);
s   = s(1);
ori = u(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function to compute the pseudo inverse. This is the same as the
% standard Matlab function, except that the default tolerance is twice as
% high.
%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision$  $Date: 2009/06/17 13:40:37 $
%   default tolerance increased by factor 2 (Robert Oostenveld, 7 Feb 2004)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = pinv(A,varargin)
[m,n] = size(A);
if n > m
  X = pinv(A',varargin{:})';
else
  [U,S,V] = svd(A,0);
  if m > 1, s = diag(S);
  elseif m == 1, s = S(1);
  else s = 0;
  end
  if nargin == 2
    tol = varargin{1};
  else
    tol = 10 * max(m,n) * max(s) * eps;
  end
  r = sum(s > tol);
  if (r == 0)
    X = zeros(size(A'),class(A));
  else
    s = diag(ones(r,1)./s(1:r));
    X = V(:,1:r)*s*U(:,1:r)';
  end
end


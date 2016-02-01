function [K,C,R,a] = prior(dims,strength,varargin)
% PRIOR used to created smoothing and shrinkage priors.
%
%   DESCRIPTION
%
%   This helper function is currently only used by dml.blogreg. See
%   dml.enet.laplacian for a related function.
%
%   input:
%   dims = the dimensions [D1 ... DN] of the measured space (e.g., X x Y x Z x T)
%   strength = n vector [I1 ... IN] specifying the strengths within each
%   dimension (the non-zero off-diagonal elements); NOTE: negative strengths give
%   positive correlations; positive strengths give negative correlations
%   and therefore more complex patterns
%
%   options:
%   'mask' is an optional boolean array specifying which index elements 
%     of the full matrix should be used as variables (default = [] = all)
%
%   'circulant' is an optional vector specifying for which dimensions we 
%     should connect the boundaries to make a torus-like shape. This has the
%     advantage that the variance is constant but forces boundary values to be the same. 
%     (default = [0 0 0 ...] = no compensation). See Gaussian Markov Random
%     field book by Rue for details.
%
%   output:
%   precision matrix K
%   covariance matrix C
%   correlation matrix R
%   average correlation between neighbouring elements a
%
%   EXAMPLE
%   [K,C,R,a] = dml.prior([28 28],[-100 -100]);
%   imagesc(reshape(R(14*28+8,:),[28 28])); axis square; colorbar
%   [K,C,R,a] = dml.prior([28 28],[100 100]);
%   imagesc(reshape(R(14*28+8,:),[28 28])); axis square; colorbar
%
%   NOTE
%   allow neighbourhood structure in the implementation
%
%   DEVELOPER
%   Marcel van Gerven (m.vangerven@donders.ru.nl)
    
  % parse options
  opt = [];
  for i=1:2:length(varargin)
    opt.(varargin{i}) = varargin{i+1};
  end
  if ~isfield(opt,'circulant'), opt.circulant = zeros(1,length(dims)); end

  nzdims = find(strength ~= 0);

  cdim = cumprod(dims);
  tdim = [1 cdim(1:(end-1))];

  if ~isfield(opt,'mask') || isempty(opt.mask)

    nfeatures = prod(dims);

    % compute the number of non-zero elements for the precision matrix
    % (diagonal plus 2 neighbours per strengthing dimension for each feature while
    % discounting pairs and ignoring boundaries plus constant term)
    K = spalloc(nfeatures,nfeatures,nfeatures+nfeatures*ceil(2^length(nzdims)/2));
    
    if any(nzdims)
      for i=1:nfeatures

        for d=nzdims

          % get neighbouring features in all dimensions

          % ignore the end boundary
          if mod(ceil(i/tdim(d)),dims(d)) ~= 0

            nbr = i + tdim(d);
            if nbr <= nfeatures
              K(i,nbr) = strength(d);
            end

          elseif opt.circulant(d)

            % compensate for the boundary
            nbr = i - (dims(d)-1)*tdim(d);
            if nbr > 0
              K(i,nbr) = strength(d);
            end

          end
        end
      end
    end
      

  else % mask specified

    opt.mask = find(opt.mask(:));
    
    nfeatures = length(opt.mask);

    nel = prod(dims);

    % compute the number of non-zero elements for the precision matrix
    % (diagonal plus 2 neighbours per dimension for each feature while
    % discounting pairs and ignoring boundaries)

    K = spalloc(nfeatures,nfeatures,ceil(1 + nfeatures+nfeatures*2^length(nzdims)/2));
    
    if any(nzdims)
      for j=1:nfeatures

        i = opt.mask(j);

        for d=nzdims

          % construct possible neighbours
          % used to detect boundaries
          pn = (i + tdim(d)):tdim(d):nel;

          % get neighbouring features in all dimensions

          % ignore the end boundary
          if any(ismember(opt.mask(opt.mask > i),pn))

            nbr = i + tdim(d);
            if nbr <= nel
              [ix,ix] = ismember(nbr,opt.mask);
              if ix
                K(j,ix) = strength(d);
              end
            end

          elseif opt.circulant(d)

            nbr = i - (dims(d)-1)*tdim(d);

            if nbr > 0
              [ix,ix] = ismember(nbr,opt.mask);
              if ix
                K(j,ix) = strength(d);
              end
            end

          end

        end
      end
    end

  end

  K = K + K';
  
  K(1:(nfeatures+1):numel(K)) = 1 - sum(K,2);
  
  if nargout > 1, C = inv(K); end
  
  if nargout > 2
    
    [m,n] = size(C);
    sigma = sqrt(diag(C)); % sqrt first to avoid under/overflow
    R = bsxfun(@rdivide,C,sigma); R = bsxfun(@rdivide,R,sigma'); % R = C ./ sigma*sigma';
    
    % Fix possible round-off problems
    t = find(abs(R) > 1); R(t) = R(t)./abs(R(t));
    R(1:m+1:end) = sign(diag(R));
    
  end
  
  if nargout > 3
    
    B = K == 0; B(1:(size(B,1)+1):end) = 0;
    a = full(mean(R(B)));

  end
  
end

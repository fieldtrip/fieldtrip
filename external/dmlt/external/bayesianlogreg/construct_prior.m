function K = construct_prior(dims,strength,varargin)
% This helper function allows the construction of the inputs for
% laplacedegenerate_ep. Bias is not yet included in this prior!
%
% input:
% dims = the dimensions [D1 ... DN] of the measured space (e.g., X x Y x Z x T)
% strength = n vector [I1 ... IN] specifying the strengths of the strengths within dimensions
%
% options:
% 'mask' is an optional vector specifying which index elements 
%    of the full matrix should be used as variables (default = [] = all)
% 'circulant' is an optional vector specifying for which dimensions we 
%    should connect the boundaries to make a torus-like shape. This has the
%    advantage that the variance of the auxiliary variables is constant but
%    forces boundary values to be the same. (default = [0 0 0 ...] = no compensation)
%

  % parse options
  opt = [];
  for i=1:2:length(varargin)
    opt.(varargin{i}) = varargin{i+1};
  end
  if ~isfield(opt,'circulant'), opt.circulant = zeros(1,length(dims)); end

  fprintf('building prior precision matrix\n');

  assert(all(strength <= 0));

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

  % make symmetric
  K = K + K';
  
  K(1:(nfeatures+1):numel(K)) = 1 - sum(K,2);
  
end

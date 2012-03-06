function KK = scale_prior(K,varargin)
% SCALE_PRIOR scales a prior precision matrix; the bias term is not yet
% part of K
%
% options:
% 'lambda' is the scaling parameter (default = 0.001)
%

  % parse options
  opt = [];
  for i=1:2:length(varargin)
    opt.(varargin{i}) = varargin{i+1};
  end
  if ~isfield(opt,'lambda'), opt.lambda = 0.001; end

  %fprintf('scaling precision matrix with lambda = %g\n', opt.lambda);

  nfeatures = size(K,1);

  K(1:(nfeatures+1):numel(K)) = 0;
  K(1:(nfeatures+1):numel(K)) = 1 - sum(K,2);
  
  [L,d1,d2] = chol(K,'lower');
   
  if d1
    error('matrix not p.d.!');
    % C = speye(nfeatures) .* inv(K);
  else    
    C = speye(nfeatures) .* fastinvc(L); % only get diagonal elements
    C     = d2*C*d2'; % reorder
  end
    
  Csqrt = sqrt(C);
  
  % KK = (Csqrt*K*Csqrt)*opt.lambda^2; % old representation in terms of regularization
  %KK = (Csqrt*K*Csqrt)/opt.lambda; % new representation in terms of scale as suggested in the paper
  KK = (Csqrt*K*Csqrt)./sqrt(opt.lambda(:)*opt.lambda(:)'); % suitable for mixed effects
  
end
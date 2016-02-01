function [beta,beta0] = elastic(X,Y,nu,lambda,options,beta,beta0)

% ELASTICNET  Elastic net implementation using coordinate descent,
%             particularly suited for sparse high-dimensional models

% X: ninput x nsamples input data
% Y: 1 x nsamples output data
% nu: ninput x 1 weights for L1 penalty
% lambda: ninput x ninput matrix for ridge penalty
% options: struct with fields
%   offset: 1 if offset is learned as well [1]
%   maxiter: maximum number of iterations [10000]
%   tol: tolerance (mean absolute difference) [1e-6]
% beta: initial beta
% beta0: initial offset
%
% beta: ninput x 1 vector
% beta0: offset ([] if not applicable)

% Parsing inputs

if nargin < 5,
    options = make_options;
else
    options = make_options(options);
end

[ninput,nsamples] = size(X);

if nargin < 4,
    lambda = 1;
end

if isscalar(lambda),
  lambda = lambda*ones(ninput,1);
end

if any(size(lambda) == 1),
  lambda = diag(lambda);
end

if nargin < 3,
    nu = 1;
end

if isscalar(nu),
    nu = nu*ones(ninput,1);
end

% Incorporating bias

nvar = ninput;
if options.offset,   
    X = [X;ones(1,nsamples)];    % expand data with 1
    nvar = nvar+1;
    lambda(nvar,nvar) = 0;
    nu = [nu;0];                 % expand nu with 0 (no regularization)
end

% Initialization

T = X*Y'/nsamples;   % nvar x 1

U = diag(lambda) + mean(X.^2,2);

if nargin < 7,
    if options.offset,
        beta0 = 0;
    else
        beta0 = [];
    end
end
if nargin < 6,
    beta = zeros(ninput,1);
end
beta = [beta;beta0];

Q = zeros(nvar,nvar);
qcomputed = zeros(nvar,1);
active = find(beta);
qcomputed(active) = 1;
Q(:,active) = X*X(active,:)'/nsamples + lambda(:,active);   % precompute Q for nonzero beta
for i=1:length(active),
    Q(active(i),active(i)) = 0;                % zero on diagonal
end
V = T-Q*beta;

% Start iterations

betaold = beta;
iter = 0;
done = 0;
while ~done,

    iter = iter+1;
    
    for i=nvar:-1:1,   % run over weights and start with offset, if applicable
        
        if abs(V(i)) <= nu(i),    % inactive
            beta(i) = 0;
        else      % active
            if V(i) > 0,
                beta(i) = (V(i) - nu(i))/U(i);
            else
                beta(i) = (V(i) + nu(i))/U(i);
            end
        end
                
        if beta(i) ~= betaold(i)
            if ~qcomputed(i),    % compute Q when needed
                Q(:,i) = X*X(i,:)'/nsamples + lambda(:,i);
                Q(i,i) = 0;
                qcomputed(i) = 1;
            end
            V = V - Q(:,i)*(beta(i) - betaold(i));
        end
        
    end
    
    
    done = (iter == options.maxiter | sum(abs(beta-betaold)) < options.tol);
    betaold = beta;
end
if iter == options.maxiter,
    fprintf('   be careful: maximum number of iterations reached\n');
end


% Parsing output

if nargout > 1,
    if options.offset,
        beta0 = beta(end);
        beta = beta(1:ninput);
    else
        beta0 = [];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%

function options = make_options(options)

if nargin < 1,
    options = struct;
end

fnames = {'offset','maxiter','tol'};
defaults = {1,1e4,1e-3};

for i=1:length(fnames),
    if ~isfield(options,fnames{i}),
        options = setfield(options,fnames{i},defaults{i});
    end
end

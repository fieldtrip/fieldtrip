function [K,W,f] = Algorithm3(Sigma, groups, lambda,normtype)
% normtype = {2,Inf}

% Groups is a vector containing the group numbers for each row
% construct cell-array with indices for faster projection
groups  = groups(:);
nGroups = max(groups);
indices = cell(nGroups,1);
for i=1:nGroups
   indices{i} = find(groups == i);
end

% ASSUMPTION: indices are contiguous!
nIndices = zeros(size(indices));
for i=1:nGroups
   nIndices(i) = length(indices{i});
end

% Get problem size
n = size(Sigma,1);

% Setup projection function
if normtype == 2
  % funProj = @projectLinf2BlockMatrix;
  funProj = @(a,b,c) projectLinf2BlockFast(a,nIndices,c); % (W, nIndices, lambda)
else
  %funProj = @projectLinf1BlockMatrix;
  funProj = @(a,b,c) projectLinf1BlockFast(a,nIndices,c); % (W, nIndices, lambda)
end

% Find initial W, using lemma 1 and diag(W) = lambda
W = initialW(Sigma,lambda(groups+nGroups*(groups-1)));
W = projectLinf1BlockMatrix(W,indices,lambda);
K = inv(Sigma + W);

% Print header
fprintf('%4s  %11s %9s %9s\n','Iter','Objective','Gap','Step');

% Main loop
i = 0; maxiter = 500; epsilon = 1e-4; alpha = 1e-3; beta = 0.5; t = 1;
f = logdet(Sigma + W,-Inf);
while (1)
   % Compute unconstrained gradient
   G = K;

   % Compute direction of step
   D = funProj(W+t*G,indices,lambda);

   % Compute initial and next objective values
f0 = f;
   ft = logdet(Sigma + D,-Inf);

   % Perform backtracking line search
   while ft < f0 + alpha * traceMatProd(D-W,G)
      % Exit with alpha is too small
      if (t < 1e-6), break; end;

      % Decrease t, recalculate direction and objective
      t = beta * t;
      D = funProj(W + t*G,indices,lambda);
      ft = logdet(Sigma + D,-Inf);
   end
   f = ft;

   % Update W and K
   W = D;
   K = inv(Sigma + W);

   % Compute duality gap
   eta = traceMatProd(Sigma,K) - n;
   for j=1:nGroups
      for l=1:nGroups

         if j==l
            y = K(indices{j},indices{l});
            mak = sum(sum(abs(y)));
         else
            y = K(indices{j},indices{l}); y = y(:);
            mak = norm(y,normtype);
         end

         if ~isempty(mak)
            eta = eta + lambda(j+nGroups*(l-1)) * mak;
         end
      end
   end

   % Increment iteration
   i = i + 1;

   % Print progress
   fprintf('%4d  %11.4e %9.2e %9.2e\n',i,f,eta,t);

   % Check stopping criterion
   if (eta < epsilon)
      fprintf('Exit: Optimal solution\n');
      break;
   elseif (i >= maxiter)
      fprintf('Exit: Maximum number of iterations reached\n');
      break;
   elseif (t < 1e-6)
      fprintf('Exit: Linesearch error\n');
      break;
   end

   % Increase t slightly
   t = t / beta;
end

end

function l = logdet(M,errorDet)

[R,p] = chol(M);
if p ~= 0
   l = errorDet;
else
   l = 2*sum(log(diag(R)));
end

global trace
if trace == 1
    global fValues
    fValues(end+1,1) = l;
    drawnow
end
end

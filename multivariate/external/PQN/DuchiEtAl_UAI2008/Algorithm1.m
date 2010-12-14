function [K,W] = Algorithm1(Sigma, lambda)

% Get problem size
n = size(Sigma,1);

% Find initial W, using lemma 1 and diag(W) = lambda
W = initialW(Sigma,diag(lambda));
K = inv(Sigma + W);

% Print header
fprintf('%4s  %11s %9s %9s\n','Iter','Objective','Gap','Step');

% Main loop
i = 0; maxiter = 1200; epsilon = 1e-4;
f = logdet(Sigma+W,-Inf);
while (1)
    % Compute unconstrained gradient
    G = K;

    % Zero components of gradient which would result in constrain violation
    G((1:n) + (0:n-1)*n) = 0; % Gii = 0
    G((W == lambda) & (G > 0)) = 0;
    G((W ==-lambda) & (G < 0)) = 0;

    % Perform line search and obtain new W
    [t,f,W] = Algorithm2(Sigma,W,K,G,lambda,f);

    % Update K
    K = inv(Sigma + W);

    % Compute duality gap
    eta = trace(Sigma * K) + sum(sum(lambda.*abs(K))) - n;

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
end
end

function [t,f,W] = Algorithm2(Sigma,W0,K,G,lambda,f0)

KG = K*G;
t  = trace(KG) / traceMatProd(KG,KG);

while(1)
    % Trial solution projected onto feasible box
    W   = W0 + t*G;
    idx = W > lambda; W(idx) =  lambda(idx); % Project - Positive part
    idx = W <-lambda; W(idx) = -lambda(idx); % Project - Negative part

    % Compute new objective
    f = logdet(Sigma + W,-Inf);
    fHat = f;
    %   D    = W - W0;
    %   KD   = K * D;
    %   fHat = f0 + trace(KD) - traceMatProd(KD,KD) / 2;

    % Test exit conditions
    if fHat >= f0, break; end;
    if t < 1e-6, break; end;

    % Reduce step length
    t = t / 2;
end

% Following is needed when using the approximate evaluation
%f = logdet(Sigma+W);

end

function l = logdet(M,errorDet)

[R,p] = chol(M);
if p ~= 0
    l = errorDet;
else
    l = 2*sum(log(diag(R)));
end;

global trace
if trace == 1
    global fValues
    fValues(end+1,1) = l;
    drawnow
end
end

function [x,f,funEvals,projects] = minConf_Sep(funObj1,funObj2,x,funProj,options)

if nargin < 5
    options = [];
end

[verbose,optTol,progTol,maxIter,memory,testOpt] = myProcessOptions(options,'verbose',2,'optTol',1e-5,...
    'progTol',1e-9,'maxIter',500,'memory',10,'testOpt',1);

% Output Log
if verbose >= 2
    if testOpt
        fprintf('%10s %10s %10s %15s %15s %15s\n','Iteration','FunEvals','Projections','Step Length','Function Val','Opt Cond');
    else
        fprintf('%10s %10s %10s %15s %15s\n','Iteration','FunEvals','Projections','Step Length','Function Val');
    end
end

% Evaluate Initial Objective
[f,g] = funObj1(x);
f = f+funObj2(x);
funEvals = 1;
projects = 0;

% Check optimality
if testOpt
optCond = max(abs(x-funProj(x-g,1)));
projects = projects+1;
if optCond < optTol
    if verbose >= 1
        fprintf('First-Order Optimality Conditions Below optTol at Initial Point\n');
    end
    return;
end
end

for i = 1:maxIter

    % Compute direction
    if i == 1
        t = min(1,1/sum(abs(g)));
        old_fvals = repmat(-inf,[memory 1]);
        old_fvals(1) = f;
        fr = f;
        alpha = 1;
    else
        y = g-g_old;
        s = x-x_old;
        alpha = (y'*s)/(y'*y);
        if alpha <= 1e-10 || alpha > 1e10
            alpha = min(1,1/sum(abs(g)));
        end
        t = 1;
        
        if i <= memory
            old_fvals(i) = f;
        else
            old_fvals = [old_fvals(2:end);f];
        end
        fr = max(old_fvals);
    end
    x_old = x;
    f_old = f;
    g_old = g;

    d = funProj(x-alpha*g,alpha)-x;
    projects = projects+1;
    
    x_new = x + t*d;
    
    [f_new,g_new] = funObj1(x_new);
    f_new = f_new + funObj2(x_new);
    funEvals = funEvals+1;

    while f_new > fr || ~isLegal(f_new)
        if verbose
        fprintf('Backtracking\n');
        end
        t = .5*t;
        
        % Check whether step has become too small
        if max(abs(t*d)) < progTol || t == 0
            if verbose == 3
                fprintf('Line Search failed\n');
            end
            t = 0;
            f_new = f;
            g_new = g;
            break;
        end
        
        x_new = x+t*d;
        [f_new,g_new] = funObj1(x_new);
        f_new = f_new + funObj2(x_new);
        funEvals = funEvals+1;
    end
    x = x_new;
    f = f_new;
    g = g_new;

    if testOpt
        optCond = max(abs(x-funProj(x-g,1)));
        projects = projects+1;
    end
       
    % Output Log
    if verbose >= 2
        if testOpt
            fprintf('%10d %10d %10d %15.5e %15.5e %15.5e\n',i,funEvals,projects,t,f,optCond);
        else
            fprintf('%10d %10d %10d %15.5e %15.5e\n',i,funEvals,projects,t,f);
        end
    end
    
    % Check Optimality
    if testOpt
    if optCond < optTol
        if verbose
            fprintf('First-order optimality below optTol\n');
        end
        break;
    end
    end
    
    if max(abs(x-x_old)) < progTol
        if verbose >= 1
            fprintf('Step size below progTol\n');
        end
        break;
    end

    if abs(f-f_old) < progTol
        if verbose >= 1
            fprintf('Function value changing by less than progTol\n');
        end
        break;
    end

    if funEvals > maxIter
        if verbose
        fprintf('Exceeded maxIter funEvals\n');
        end
        break
    end
end
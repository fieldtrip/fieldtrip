function [w] = L1groupGraft(funObj,w,groups,lambda,options)

if nargin < 5
    options = [];
end

[maxIter,optTol] = myProcessOptions(options,'maxIter',500,'optTol',1e-6);

nVars = length(w);
nGroups = max(groups);

reg = sqrt(accumarray(groups(groups~=0),w(groups~=0).^2));

% Compute Initial Free Variables
free = ones(nVars,1);
for s = 1:nGroups
    if reg(s) == 0
        free(groups==s) = 0;
    end
end

% Optimize Initial Free Variables
subOptions = options;
subOptions.TolFun = 1e-16;
subOptions.TolX = 1e-16;
subOptions.Display = 'none';
if any(free == 1)
    [w(free==1) f junk2 output] = minFunc(@subGradient,w(free==1),subOptions,free,funObj,groups,lambda);
    fEvals = output.funcCount;
end

i = 1;
maxViolGroup = -2;
old_maxViolGroup = -1;
while fEvals < maxIter
    reg = sqrt(accumarray(groups(groups~=0),w(groups~=0).^2));
    [f,g] = subGradient(w,ones(nVars,1),funObj,groups,lambda);

    fprintf('%5d %5d %15.5f %15.5e %5d\n',i,fEvals,f,sum(abs(g)),sum(reg > 0));

    if sum(abs(g)) < optTol
        fprintf('Solution Found\n');
        break;
    end
    
    if all(free==1)
        fprintf('Stuck\n');
    end
    
    % Compute Free Variables
    free = ones(nVars,1);
    for s = 1:nGroups
        if reg(s) == 0
            free(groups==s) = 0;
        end
    end

    % Add Group with biggest sub-gradient
    gradReg = sqrt(accumarray(groups(groups~=0),g(groups~=0).^2));
    [maxViol maxViolGroup] = max(gradReg);
    free(groups==maxViolGroup) = 1;
    
    if maxViolGroup == old_maxViolGroup
        fprintf('Stuck (optTol = %15.5e)\n',sum(abs(g)));
        break;
    end
    old_maxViolGroup = maxViolGroup;

    % Optimize Free Variables
    if any(free == 1)
        [w(free==1) f junk2 output] = minFunc(@subGradient,w(free==1),subOptions,free,funObj,groups,lambda);
        fEvals = fEvals+output.funcCount;
    end

    i = i + 1;
end

end

function [f,g] = subGradient(wSub,free,funObj,groups,lambda)
w = zeros(size(free));
w(free==1) = wSub;

[f,g] = funObj(w);
f = f + lambda*sum(sqrt(accumarray(groups(groups~=0),w(groups~=0).^2)));
g = getGroupSlope(w,lambda,g,groups,1e-4);
g = g(free==1);
end
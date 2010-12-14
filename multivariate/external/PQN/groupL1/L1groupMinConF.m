function [w,f] = L1groupMinConF(funObj,w,groups,lambda,options)
% [w] = L1groupMinConF(funObj,w,groups,lambda,options)

if nargin < 5
    options = [];
end

[normType,mode,optTol] = myProcessOptions(options,'normType',2,'mode','spg','optTol',1e-6);

nVars = length(w);
nGroups = max(groups);

% Make initial values for auxiliary variables
wAlpha = [w;zeros(nGroups,1)];
for g = 1:nGroups
    if normType == 2
        wAlpha(nVars+g) = norm(w(groups==g));
    else
        if any(groups==g)
            wAlpha(nVars+g) = max(abs(w(groups==g)));
        end
    end
end

% Make Objective Function
wrapFunObj = @(w)auxGroupLoss(w,groups,lambda,funObj);

if normType == 2
    switch mode
        case 'barrier'
            wAlpha(nVars+1:end) = wAlpha(nVars+1:end)+1e-2;
            funCon = @(w)groupL2Barrier(w,groups);
            [wAlpha,f,fEvals] = minConF_LB(wrapFunObj,wAlpha,funCon,options);
        case 'penalty'
            % Doesn't work
            funCon = @(w)groupL2Penalty(w,groups);
            [wAlpha,f,fEvals] = minConF_Pen(wrapFunObj,wAlpha,funCon,options);
        case 'spg'
            [groupStart,groupPtr] = groupl1_makeGroupPointers(groups);
            funProj = @(w)auxGroupL2Project(w,nVars,groupStart,groupPtr);
            %funProj = @(w)auxGroupL2Proj(w,groups);
            wAlpha = minConF_SPG(wrapFunObj,wAlpha,funProj,options);
        case 'sop'
            [groupStart,groupPtr] = groupl1_makeGroupPointers(groups);
            funProj = @(w)auxGroupL2Project(w,nVars,groupStart,groupPtr);
            options.bbInit = 0;
            options.SPGoptTol = 1e-6;
            options.SPGiters = 500;
            options.maxProject = inf;
            options.SPGtestOpt = 1;
            [wAlpha,f] = minConF_PQN(wrapFunObj,wAlpha,funProj,options);
        case 'interior'
            wAlpha(nVars+1:end) = wAlpha(nVars+1:end)+1e-2;
            funCon = @(w,lambda)groupL2Residuals(w,lambda,groups);
            wAlpha = minConF_IP2(wrapFunObj,wAlpha,funCon,options);
        case 'sep'
            funObj1 = @(w)funObj(w);
            funObj2 = @(w)groupL1regularizer(w,lambda,groups);
            funProj = @(w,stepSize)groupSoftThreshold(w,stepSize,lambda,groups);
            wAlpha = minConF_Sep(funObj1,funObj2,w,funProj,options);
    end
else
   switch mode 
       case 'barrier'
           [A,b] = makeL1LinfConstraints(groups);
           wAlpha(nVars+1:end) = wAlpha(nVars+1:end)+1e-2;
           funCon = @(w)linearBarrier(w,A,b);
           [wAlpha,f,fEvals] = minConF_LB(wrapFunObj,wAlpha,funCon,options);
       case 'penalty'
           [A,b] = makeL1LinfConstraints(groups);
           funCon = @(w)linearInequalityPenalty(w,A,b);
           [wAlpha,f,fEvals] = minConF_Pen(wrapFunObj,wAlpha,funCon,options);
       case 'spg'
           funProj = @(w)auxGroupLinfProj(w,groups);
           wAlpha = minConF_SPG(wrapFunObj,wAlpha,funProj,options);
       case 'sop'
           funProj = @(w)auxGroupLinfProj(w,groups);
           funSOP = @(w,g,H)SOP_SPG(w,g,H,funProj);
           wAlpha = minConF_SOP(wrapFunObj,wAlpha,funProj,funSOP,options);
       case 'active'
           % Doesn't work (constraints are degenerate)
           [A,b] = makeL1LinfConstraints(groups);
           wAlpha = minConF_AS(wrapFunObj,wAlpha,A,b,options);
       case 'interior'
           [A,b] = makeL1LinfConstraints(groups);
           wAlpha(nVars+1:end) = wAlpha(nVars+1:end)+1e-2;
           wAlpha = minConF_IP(wrapFunObj,wAlpha,-A,b,options);
   end
end
w = wAlpha(1:nVars);

end

function [A,b] = makeL1LinfConstraints(groups)
nVars = length(groups);
nGroups = max(groups);
A = zeros(0,nVars+nGroups);
j = 1;
for i = 1:nVars
    if groups(i) ~= 0
        A(j,i) = 1;
        A(j,nVars+groups(i)) = 1;
        A(j+1,i) = -1;
        A(j+1,nVars+groups(i)) = 1;
        j = j+2;
    end
end
b = zeros(size(A,1),1);
end

function [f] = groupL1regularizer(w,lambda,groups)
f = lambda*sum(sqrt(accumarray(groups(groups~=0),w(groups~=0).^2)));
end

function [w] = groupSoftThreshold(w,alpha,lambda,groups)
    nGroups = max(groups);
    reg = sqrt(accumarray(groups(groups~=0),w(groups~=0).^2));
    for g = 1:nGroups
        w(groups==g) = (w(groups==g)/reg(g))*max(0,reg(g)-lambda*alpha);
    end
end
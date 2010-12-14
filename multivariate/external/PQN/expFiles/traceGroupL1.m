function [f,g] = traceGroupL1(w,lambda,groups,funObj)
lambdaL2 = 1e-4;

[f,g] = funObj(w);
f = f + (lambdaL2/2)*sum(w.^2);
g = g + lambdaL2*w;

global fValues;

fValues(end+1,1) = f + lambda*sum(sqrt(accumarray(groups(groups~=0),w(groups~=0).^2)));

drawnow
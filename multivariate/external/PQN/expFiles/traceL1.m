function [f,g] = traceL1(w,lambda,funObj)
[f,g] = funObj(w);

global fValues;
fValues(end+1,1) = f + sum(lambda.*abs(w));

drawnow
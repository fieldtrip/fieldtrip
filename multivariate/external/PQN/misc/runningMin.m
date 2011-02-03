function [minVals] = runningMin(x,y)
% Compute running min of y, and use these indexes to choose elements of x
% (typically x is the same y)

if nargin < 2
    y = x;
end

minVals = zeros(length(x),1);
for i = 1:length(x)
    [junk minPos] = min(y(1:i));
    minVals(i) = min(x(minPos));
end
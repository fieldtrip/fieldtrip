function [retval] = rv(d1, d2);

% RV returns the relative residual variance between measured and simulated data
%
% rv = rv(measured, simulated)

% Copyright (C) 1999, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

retval = sum((d1-d2).^2) ./ sum(d1.^2);


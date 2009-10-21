function [retval] = rv(d1, d2);

% RV returns the relative residual variance between measured and simulated data
%
% rv = rv(measured, simulated)

% Copyright (C) 1999, Robert Oostenveld
%
% $Log: rv.m,v $
% Revision 1.2  2003/03/11 14:45:37  roberto
% updated help and copyrights
%

retval = sum((d1-d2).^2) ./ sum(d1.^2);


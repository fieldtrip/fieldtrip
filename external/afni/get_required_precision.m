function p=get_required_precision(d,epsilon)
% Finds the number of decimal positions required to print all data
% elements.
%
% P=GET_REQUIRED_PRECISION(D) returns the least number of decimals
% 'after the period' to represent all numbers D without loss of precision.
% D should be an array with real numbers and can be of any size.
%
% P=GET_REQUIRED_PRECISION(D,EPSILON) allows for an imprecision of
% at most EPSILON.
%
% In the current implementation, this function always returns a value
% between 0 and 7 (inclusive).
%
% Examples: GET_REQUIRED_PRECISION([2 5])       returns 0
%           GET_REQUIRED_PRECISION([2 5.532])   returns 3
%
% NNO Jan 2010

if nargin<2
    epsilon=0; % we want high precision of course
end

% max precision is 7, because that seems more or less the limit of 32 bit
% floats

for p=0:7
    if max(abs(rem(d(:),10^(-p))))<=epsilon
        return
    end
end
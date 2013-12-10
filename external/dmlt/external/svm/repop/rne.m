function [Z]=rne(X,Y)
% 
% Replicating not equal
%
% Does element by element operations on X and Y where non-same sized
% dimensions are implicity wrapped round to match the size of the larger
% to give a result matrix Z with size max(size(X),size(Y));
%
% In this case returns logical array with true where X~=Y and 0 otherwise.
%
% N.B. for complex inputs this compares the *norms* of the values, *not*
%      the real parts as GE does!
% See also repops, ne, eq
%
% Copyright 2006- by Jason D.R. Farquhar (jdrf@zepler.org)
% Inspired by code by Douglas M. Schwarz & Aki Vehtari.
Z=repop(X,'~=',Y);

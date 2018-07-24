function K=linKer(X,Z,dir)
% Compute linear kernel using tprod
%
% K=linKer(X,Z,dir)
%undefined reference to `copymxInfo'
% Inputs:
% X   -- [N1 x ... ] or [ ... x N1 ] data matrix
% Y   -- [N2 x ... ] or [ ... x N2 ] data matrix
% dir -- flag to indicate the trials dimension.
%        1  -> trials in the first dimension, i.e.X= [N x ... ]
%        -1 -> trials in the last dimension, i.e. X= [... x N ]
% Outputs
% K   -- [N1 x N2] linear kernel matrix
%
% Copyright 2006-     by Jason D.R. Farquhar (jdrf@zepler.org)

% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied

if ( nargin < 2 ) Z=X; end;
if ( nargin < 3 ) dir=1; end;
if ( dir==-1 ) && (ndims(X)>2 || ndims(Z)>2)% dir indicates if row or col vectors
   K=tprod(X,[-(1:ndims(X)-1) 1],Z,[-(1:ndims(Z)-1) 2],'n');
elseif ( dir==-1)
   K=transpose(X)*Z;
elseif ( dir==1 ) && (ndims(X)>2 || ndims(Z)>2)
   K=tprod(X,[1 -(2:ndims(X))],Z,[2 -(2:ndims(Z))],'n');
elseif ( dir==1 )
   K=X*transpose(Z);
end
function varargout = spm_krutil(varargin)
% Men doing things together (with Kronecker Tensor Products)
%
% This is a multi-format function where the number of input and output
% arguments determine which precise function is being accessed.  It's
% general purpose is to speed up certain calculations based on
% Kronecker Tensor products. 
%
% Some of the functions are based on the fact that kron(A'*A,B'*B) 
% is DRAMATICALLY faster to calculate than kron(A,B)'*kron(A,B).
% See Ashburner & Friston, 1999, Human Brain Map 7:254-266.
%
% Other functions are simply slight improvements on the details of the
% Matlab implementation of kron such that they are a bit faster and
% less memory hungry.
%
%_______________________________________________________________________
% Format of cunning functions of first kind.
%
% FORMAT beta = spm_krutil(img,Bx,By,0)
%
% Equivalent to: beta = kron(By,Bx)'*img(:)
%
% FORMAT alpha = spm_krutil(img,Bx,By,1)
%
% Equivalent to: alpha = kron(By,Bx)'*diag(img(:))*kron(By,Bx)
%
% FORMAT alpha = spm_krutil(img,B1x,B1y,B2x,B2y)
%
% Equivalent to: alpha = kron(B1y,B1x)'*diag(img(:))*kron(B2y,B2x)
%
%_______________________________________________________________________
% Format of functions of the second kind.
%
% FORMAT C = spm_krutil(A,B)
%
% is equivalent to, but faster than, C = kron(A,B)
%
% FORMAT spm_krutil(A,B,C)
%
% is equivalent to, but faster and less memory consuming than
% C = C + kron(A,B). N.B. this latter form has the slightly dangerous
% property of changing one of its input arguments. Be very careful
% when using it like this.
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner & Jesper Andersson
% $Id: spm_krutil.m 1143 2008-02-07 19:33:33Z spm $


%-This is merely the help file for the compiled routine
error('spm_krutil.c not compiled - see Makefile')

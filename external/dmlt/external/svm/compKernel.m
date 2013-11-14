function K = compKernel(X,Z,kerType,varargin)
%K = compKernel(X,Z,kerType,varargin)
%
%Inputs:
% X: Input objects: [N1 x d]
% Z: Input objects: [N2 x d]
% kerType: kernel type. this is either an actual kernel matrix, or String.
%      'linear'  - Linear      K(i,j) = X(i,:)*Z(j,:)
%      'poly'    - Polynomial  K(i,j) = (X(i,:)*Z(j,:)+par(2))^par(1)
%      'rbf'     - RBF         K(i,j) = exp(-|X(i,:)-Z(j,:)|^2/(2*par))
%      'precomp' - par holds kernel            K(i,j) = par(i,j)
%      @kerfn    - function_handle to function 'func' such that :
%                      K(i,j)=func(X,Z,varargin)
%
% N.B. use nlinear, npoly, nrbf, etc. to compute the normalised kernel,
%      i.e. which has all ones on the diagonal.
%
%   varargin: parameter(s) of the kernel.
%
%Outputs:
% K: the computed kernel matrix.
%
% Modified : Jason Farquhar (jdrf@zepler.org) 2005-08-01
% Version:  1.0.1
%
% Copyright 2006-     by Jason D.R. Farquhar (jdrf@zepler.org)

% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied


if (nargin < 3) % check correct number of arguments
  error('Insufficient arguments'); return;
end
if ( isinteger(X) ) X=single(X); end % need floats for products!
% empty Z means use X
isgram=false;
dir   =1;  % Assume row vectors in X
if ( isempty(Z) ) Z=X; isgram=true; elseif ( isinteger(Z) ) Z=single(Z); end;

if ( isstr(kerType) )
 switch lower(kerType)

  case {'linear','nlinear'}; % linear
   K = linKer(X,Z,dir);

  case {'poly','npoly'};     % polynomial
   if(numel(varargin)<2) varargin{2}=1;end;
   K = (linKer(X,Z,dir)+varargin{2}).^varargin{1};

  case {'rbf','nrbf'};       % Radial basis function, a.k.a. gaussian
   K = sqDist(X,Z,dir); % pairwise distance
   K = exp( - .5 * K / varargin{1} ) ;
   if(isgram) K=.5*(K+K'); end % to avoid rounding error problems....

  case 'precomp';            % Given in the input
   if ( all( size(varargin) == [size(X1,1) size(X2,1)] ) )
     K=varargin{1};
   else
     error('Kernel matrix does not match input dimensions');
   end

  otherwise;
   if ( exist(kerType)>1 )   % String which specifies function on the path
      K=feval(kerType,X,Z,varargin{:});
   else
     error(['Unrecognised kernel type : ' kerType]);
   end
 end

elseif ( isa(kerType,'function_handle') ) % function handle
 % use this handle
 K=feval(kerType,X,Z,varargin{:});
else
  error('Unknown kernel type');return
end

if ( isequal(kerType(1),'n') ) % normalise computed kernel
  if ( isgram ) % Need to compute K(X,X)_i,i and K(Z,Z)_j,j
     Nr = diag(K); Nc = diag(K);
  else
     for i=1:size(X,1); Nr(i,1)=compKernel(X(i,:),[],kerType,varargin{:});end;
     for i=1:size(Z,1); Nc(i,1)=compKernel(Z(i,:),[],kerType,varargin{:});end;
  end
  K = repop(repop(K,sqrt(Nr),'./'),sqrt(Nc)','./');
end;

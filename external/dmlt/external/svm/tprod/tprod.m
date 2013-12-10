function [varargout]=tprod(varargin)
% Multi-dimensional generalisation of matrix multiplication
%
% function [Z]=tprod(X,xdimspec,Y,ydimspec[,optStr,blksz])
%
% This function computes a generalised multi-dimensional matrix product
% based upon the Einstein Summation Convention (ESC).  This means
% given 2 n-d inputs:
%   X = [ A x B x C x E .... ], denoted in ESC as X_{abce}
%   Y = [ D x F x G x H .... ], denoted in ESC as Y_{dfgh}
% we define a particular tensor operation producing the result Z as
%  Z_{cedf} = X_{abce}Y_{dfab}
% where, we form an inner-product (multiply+sum) for the dimensions with
% *matching* labels on the Right Hand Side (RHS) and an outer-product over
% the remaining dimensions.  Note, that in conventional ESC the *order* of
% dimensions in Z is specified on the LHS where the matching label specifies
% its location.
% 
% Tprod calls closely follow this convention, the tprod equivalent of the
% above is[1]:
%  Z = tprod(X,[-1 -2 1 2],Y,[3 4 -1 -2])
% here, *matching negatively* labelled RHS dimensions are inner-product
% dimensionss, whilst the *positively* labelled dimensions directly
% specify the position of that dimension in the result Z. Hence only 2
% dimension-specifications are needed to unambigiously specify this
% tensor product.  
%
% [1] Note: if you find this syntax to difficult the ETPROD wrapper function
% is provided to directly make calls in ESC syntax, e.g.
%    Z = etprod('cedf',X,'abce',Y,'dfab');
%
% It is perhaps easiest to understand the calling syntax with some simple
% examples: 
% Simple 1-d cases
%   X = randn(100,1);                      % make 100 1-d data points
%   Z = tprod(X,1:ndims(X),X,1:ndims(X));           % element-wise multiply
%   Z = tprod(X,1:ndims(X),X,ndims(X)+1:ndims(X));  % outer-product
%   Z = tprod(X,-1,X,-1);                           % inner product   
% Some 2-d statistical examples
%   X=randn(10,100);                       % 10d points x 100 trials
%   Z=tprod(X,[-1 2],X,[-1 2]);            % squared norm of each trial
%   Z=tprod(X,[1 -2],X,[2 -2]);            % covariance over trials
% More complex 3-d cases
%   X = randn(10,20,100);          % dim x samples x trials set of timeseries
%   sf= randn(10,1);               % dim x 1 spatial filter
%   tf= randn(20,1);               % samples x 1 temporal filter
%   Z=tprod(X,[-1 1 2],sf,[-1 1]); % spatially filter Z -> [1 x samp x trials ]
%   Z=tprod(X,[1 -2 3],tf,[-2 2]); % temporaly fitler Z -> [dim x 1 x trials ]
%   % OP over dim, IP over samples, sequ over trial = per trial covariance
%   Z=tprod(X,[1 -2 3],X,[2 -2 3])/size(X,3); 
%
% INPUTS:
%  X        - n-d double/single matrix
%
%  xdimspec - signed label for each X dimension. Interperted as:
%              1) 0 labels must come from singlenton dims and means they are
%              squeezed out of the input.
%              2) NEGATIVE labels must come in matched pairs in both X and Y 
%                 and denote inner-product dimensions
%              3) POSITIVE labels denote the position of this dimension in
%              the output matrix Z.  Positive labels must be unique in X.
%              Depending on whether the same label occurs in Y 2 conditions
%              can occur:
%                a) X label has NO match in Y.  Then this is an
%                outer-product dimension.  
%                b) X label matches a label in Y. Then this dimension is an
%                aligned in both X and Y, such that they increment together
%                -- as if there was an outer loop over these dims indicies
%
%  Y        - m-d double/single matrix, 
%             N.B. if Y==[], then it is assumed to be a copy of X.
%
%  ydimspec - signed label for each Y dimension.  If not given yaccdim
%             defaults to -(1:# negative labels in (xdimspec)) followed by
%             enough positive lables to put the remaining dims after the X
%             dims in the output. (so it accumlates the first dims and
%             outer-prods the rest)
%
%  optStr  - String of single character control options,
%            'm'= don't use the, fast but perhaps not memory efficient,
%                 MATLAB code when possible. 
%            'n'=use the *new* calling convention: tprod(X,xdimspec,Y,ydimspec)
%            'o'=use the *old* calling convention: tprod(X,Y,xdimspec,ydimspec)
%  blksz    - Internally tprod computes the results in blocks of blksz size in 
%             order to keep information efficiently in the cache.  TPROD 
%             defaults this size to a size of 16, i.e. 16x16 blocks of doubles.
%             On different machines (with different cache sizes) tweaking this
%             parameter may result some speedups.
% 
% OUTPUT:
%  Z       - n-d double matrix with the size given by the sizes of the
%            POSITIVE labels in xdimspec/ydimspec
%
% See Also:  tprod_testcases, etprod
%
% Class support of input:
%     float: double, single
% 
% 
% Copyright 2006-     by Jason D.R. Farquhar (jdrf@zepler.org)

% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied


% The rest of this code is a mex-hiding mechanism which compilies the mex if
% this runs and recursivly calls itself.  
% Based upon code from: http://theoval.sys.uea.ac.uk/matlab
cwd  = pwd; % store the current working directory
name = mfilename('fullpath'); % get the directory where we are
% find out what directory it is defined in
name(name=='\')='/';
dir=name(1:max(find(name == '/')-1));
try % try changing to that directory
   cd(dir);
catch   % this should never happen, but just in case!
   cd(cwd);
   error(['unable to locate directory containing ''' name '.m''']);
end

try % try recompiling the MEX file
   fprintf(['Compiling ' mfilename ' for first use\n']);
   mex('ddtprod.c','dstprod.c','sdtprod.c','sstprod.c','tprod_util.c','tprod_mex.c','mxInfo.c','mxInfo_mex.c','-O','-output',mfilename);
   fprintf('done\n');
catch
   % this may well happen happen, get back to current working directory!
   cd(cwd);
   error('unable to compile MEX version of ''%s''%s\n%s%s', name, ...
         ', please make sure your', 'MEX compiler is set up correctly', ...
         ' (try ''mex -setup'').');
end

cd(cwd); % change back to the current working directory
rehash;  % refresh the function and file system caches

% recursively invoke MEX version using the same input and output arguments
[varargout{1:nargout}] = feval(mfilename, varargin{:});

% bye bye...

return;

function [varargout]=tprod(varargin)
% function [Z]=tprod(X,Y,xdimspec,ydimspec[,optStr,blksz])
%
% Multi-dimensional generalisation of matrix multiplication
%
% This function computes a generalised multi-dimensional matrix product
% based upon the einstein summation convention (plus extras).  This means
% given 2 n-d inputs:
%   X = [ A x B x C x E .... ] 
%   Y = [ D x F x G x H .... ]
% we define the result, Z, to be
%  Z = X_{-1,-2,1,2}Y_{3,4,-1,-2} = [ C x E x D x F ]
%
% N.B. if Y==[], then it is assumed to be a copy of X.
% This result is produced by forming an inner-product (multiply+sum) for the
% pairs of dimensions which have the same (negative) label (i.e. -1 => X dim
% 1 and Y dim 3, and -2 => X dim 2 and Y dim 4) and forming an outer product
% for all the remaining dimensions, with the positive label determing where
% this dimension is placed in the output.
%
% Depending on the state of the squeeze flag the size of Z is given by the 
% size of X followed by the size of Y but:
%   with the accumulated dimensions removed with squeeze turned on (default)
%    or, 
%   with the accumulated dimensions size set to 1 when squeeze is turned off.
%
% Examples:
%   X = randn(100,1);       % 100d data point
%   Z = tprod(X,[],1:ndims(X),1:ndims(X));        % outer product
%   Z = tprod(X,[],-1,-1);      % inner product   
%   X = randn(100,100,100); % dim x samp x trials set of timeseries
%   sf= randn(100,1);       % dim x 1 spatial filter
%   tf= randn(100,1);       % dim x 1 temporal filter
%   Z=tprod(X,sf,[-1 1 2],[-1 1]); % spatially filter Z -> [ 1 x samp x trials ]
%   Z=tprod(X,tf,[1 -2 3],[-2 2]); % temporaly fitler Z -> [ dim x 1 x trials ]
%
%   X=randn(10,100);           % 10d x 100 trial
%   Z=tprod(X,X,[-1 2],[1 2]); % sum squared value over trial
%   Z=tprod(X,X,[0 -2],[0 2]); % covariance for each trial
%   X=randn(10,100,400);       % dim x samp x trials set of timeseries
%   % OP over dim, IP over samp, sequ over trial = per trial dim covariance
%   Z=tprod(X,X,[0 2 -3],[0 2 3])/size(X,3); 
%
% INPUTS:
%  X        - n-d double matrix
%  Y        - n-d double matrix
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
%                with an inner tprod over the remaining non-aligned
%                dimensions.
%  ydimspec - signed label for each Y dimension.  If not given yaccdim
%             defaults to -(1:# negative labels in (xdimspec)) followed by
%             enough positive lables to put the remaining dims after the X
%             dims in the output. (so it accumlates the first dims and
%             outer-prods the rest)
%  optStr  - String of options to turn *OFF* functions, can contain be either:
%            'm'= don't use the, fast but perhaps not memory efficient,
%                 MATLAB code when possible. 
%  blksz    - Internally tprod computes the results in blocks of blksz size in 
%             order to keep information efficiently in the cache.  TPROD 
%             defaults this size to a size of 16, i.e. 16x16 blocks of doubles.
%             On different machines (with different cache sizes) tweaking this
%             parameter may result some significant speedups.
% 
%
% OUTPUT:
%  Z       - n-d double matrix with the size given by the sizes of the
%            POSITIVE labels in xdimspec/ydimspec
%
% See Also:  tprod_testcases, etprod
%
% Class support of input P:
%     float: double
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
   mex('ddtprod.c','dstprod.c','sdtprod.c','sstprod.c','tprod_util.c','tprod_mex.c','mxInfo.c','-O','-output',mfilename);
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

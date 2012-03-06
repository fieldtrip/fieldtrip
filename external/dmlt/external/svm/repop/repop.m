function [varargout]=repop(varargin)
% Replicating arithemetical and logical operators.
%
% [Z]=repop(X,operator,Y [,options])
%
% Does element by element operations on X and Y where non-same sized
% dimensions are implicity wrapped round to match the size of the larger
% to give a result matrix Z with size max(size(X),size(Y));
%
% What this means is that if you give a [Nx1] and a [1xM] input you get 
% a [NxM] output, or if you give it a [Nx1xQ] and a [NxMx1xP] you get 
% [NxMxQxP] output etc..  
% Note that the size of the input *completely and uniquely* determines the 
% size of the output.
%
% In general this is at least 2x faster than the equivalent matlab code
% using repmats and has the advantage of requiring no additional memory.
%
% Example Usage:
%     X = randn(10000,10);                  % example signal with data in rows
%     stdX = repop(X,'-',mean(X,1));        % subtract the mean vector
%     stdX = repop(stdX,'/',std(stdX,0,1)); % divide by std-deviation
%
% Operator can be one of:
%
% Arthemetical -- returns a double matrix
%   '+','.+',plus   - Implicitly repmatted elementwise addition
%   '-','.-',minus  - Implicitly repmatted elementwise addition
%   '*','.*',times  - Implicitly repmatted elementwise multiplication
%   '^','.^',power  - Implicitly repmatted elementwise raise X to power Y
%   '\','.\',ldivide- Implicitly repmatted elementwise divide Y by X
%   '/','./',rdivide- Implicitly repmatted elementwise divide X by Y
%   'min'           - Implicitly repmatted elementwise min of X by Y
%   'max'           - Implicitly repmatted elementwise max of X by Y
%
% Relational -- returns a logical matrix
% N.B. for complex inputs the <,>,<=,>= operators are based upon abs(x) 
% (not real(x) as in matlab)
%   '==',eq         - Implicitly repmatted elementwise equality
%   '~=',ne         - Implicitly repmatted elementwise dis-equality
%   '<' ,lt         - Implicitly repmatted elementwise less than
%   '>' ,gt         - Implicitly repmatted elementwise greater than
%   '<=',le         - Implicitly repmatted elementwise less than equal
%   '>=',ge         - Implicitly repmatted elementwise greater than equal
%
% N.B. the operator can go in any of the 3 argument positions, i.e.
% these are all valid: repop(X,'-',Y), repop(X,Y,'-'), repop('-',X,Y)
%
% The optional final argument is a string of single letter switches
% consisting off
%  'm'  -- allows replication of non-unit dimensions if the larger
%          dimensions size is an integer multiple of the smaller ones
%  'n'  -- allow replication of non-unit dimensions in *all* cases
%  'i'  -- perform "inplace" operation.  This means that we use a
%          *dangerous* matlab *hack* to perform the operation *without*
%          allocating new memory for the output, but by simply overwriting
%          the memory used for X in the input.  Thus the following code is
%          a memory (and time) efficient way to increment X.
%              X = repop(X,'+',1,'i');
%
%
% Class support of input P:
%     float: double, single
% 
% SEE ALSO:   repop_testcases rplus rminus rtimes rpower rldivide rrdivide
%             req rne rlt rgt rle rge
%
% Copyright 2006-     by Jason D.R. Farquhar (jdrf@zepler.org)
% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied
%
% Inspired by code from Douglas M. Schwarz and & Aki Vehtari.

% The rest of this code is a mex-hiding mechanism which compilies the mex if
% this runs and recursivly calls itself.  
% Based upon code from: http://theoval.sys.uea.ac.uk/matlab
cwd  = pwd; % store the current working directory
name = mfilename('fullpath'); % get the directory where we are
% find out what directory it is defined in
name(name=='\')='/'; % deal with dos'isms
dir=name(1:max(find(name == '/')-1)); % dir is everything before final '/'
try % try changing to that directory
   cd(dir);
catch   % this should never happen, but just in case!
   cd(cwd);
   error(['unable to locate directory containing ''' name '.m''']);
end

try % try recompiling the MEX file
   fprintf(['Compiling ' mfilename ' for first use\n']);
   mex('ddrepop.c','dsrepop.c','sdrepop.c','ssrepop.c','repop_util.c','repop_mex.c','mxInfo.c','mxInfo_mex.c','-O','-output',mfilename);
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

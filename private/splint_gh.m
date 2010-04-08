function [varargout] = funname(varargin)

% SPLINT_GH implements equations (3) and (5b) of Perrin 1989
% for simultaneous computation of multiple values

% Copyright (C) 2004-2009, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

% compile the missing mex file on the fly
% remember the original working directory
pwdir = pwd;

% determine the name and full path of this function
funname = mfilename('fullpath');
mexsrc  = [funname '.c'];
[mexdir, mexname] = fileparts(funname);

try
  % try to compile the mex file on the fly
  warning('trying to compile MEX file from %s', mexsrc);
  cd(mexdir);
  mex(mexsrc);
  cd(pwdir);
  success = true;

catch
  % compilation failed
  disp(lasterr);
  error('could not locate MEX file for %s', mexname);
  cd(pwdir);
  success = false;
end

if success
  % execute the mex file that was juist created
  funname   = mfilename;
  funhandle = str2func(funname);
  [varargout{1:nargout}] = funhandle(varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE FOLLOWING CODE CORRESPONDS WITH THE ORIGINAL IMPLEMENTATION
% function [gx, hx] = gh(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% M = 4;          % constant in denominator
% N = 9;          % number of terms for series expansion
% p  = zeros(1,N);
% gx = zeros(size(x));
% hx = zeros(size(x));
% x(find(x>1)) = 1;       % to avoid rounding off errors
% x(find(x<-1)) = -1;     % to avoid rounding off errors
% for i=1:size(x,1)
% for j=1:size(x,2)
%   for k=1:N
%     p(k) = plgndr(k,0,x(i,j));
%   end
%   gx(i,j) =  sum((2*(1:N)+1) ./ ((1:N).*((1:N)+1)).^M     .* p) / (4*pi);
%   hx(i,j) = -sum((2*(1:N)+1) ./ ((1:N).*((1:N)+1)).^(M-1) .* p) / (4*pi);
% end
% end
%

function f = str2fun(f)
% STR2FUN - Compatibility wrapper to str2func
%
%           Description
%           the code of function:
%
%           if exist('str2func')
%             f=str2func(f);
%           end

% Copyright (c) 1998-2004 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if exist('str2func')
  f=str2func(f);
end

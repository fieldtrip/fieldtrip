function [x,rec] = fsearch(fun,x0,opt,varargin)
% FSEARCH Finds the minimum of a combinatorial function using forward search
%
%   X = FSEARCH(FUN, X0) attempts to find a combination of elements
%   of X0 which locally minimize the function FUN using forward
%   search. FUN accepts input X and returns scalar function value F
%   evaluated at X. X0 must be a vector of candidate indexes. 
%   Returned X contains indexes locally minimizing the function
%   FUN.
%
%   X = FSEARCH(FUN, X0, OPTIONS) allows use of optional
%   search parameters. See FSEARCH_OPT for details.
%
%   X = FSEARCH(FUN, X0, OPTIONS, P1, ..., Pn) P1,...,Pn are
%   additional parameters passed to the function FUN.
%
%   [X,REC] = FSEARCH(FUN, X0, ...) returns record of search as a
%   array of structs.
%
%   In forward search the elements are chosen one at the time. The
%   element that minimizes the function value when added to the set
%   is chosen at each round. In a case of tie choose randomly.
%

% Copyright (c) Aki Vehtari (2007)

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% 2007-02-23  Aki Vehtari  <Aki.Vehtari@hut.fi>
%             Use same parameter order as in fminunc etc.
%             Minimize instead of maximize
%             Total rewrite.

nx = size(x0,2);  % number of elements
chosen = [];      % the elements chosen sofar
nchosen = 0;      % number of chosen elements
% Options
opt = fsearch_opt(opt);

% The loop     
value=Inf;
minvalue=value;
while (nchosen < opt.nsel)
  values=zeros(nx,1);
  for i1=1:nx
    values(i1) = feval(fun,union(chosen,x0(i1)),varargin{:});
    if (opt.display >= 2)
      fprintf('%d %.5g\n', x0(i1),  values(i1));
    end
  end
  % Let's get the index of the best element
  % In a case of tie, choose randomly
  value=min(values);
  mini=randpick(find(values==value));
  if (opt.stop && value>=minvalue)
    if opt.display
      fprintf('Stopping because value not decreasing\n');
    end
    break
  end
  chosen = union(chosen,x0(mini));   % Add the element into the set of chosen
  nchosen = nchosen + 1;
  % copy some information of this round to struct rec
  rec(nchosen).chosen = chosen;
  rec(nchosen).candidates = x0;
  rec(nchosen).values = values;
  x0(mini) = [];                     % remove it from the x0
  nx=nx-1;
  if value<minvalue
    minvalue=value;
    minvaluei=nchosen;
  end

  if (opt.display >= 1)
    fprintf(' Value: %.4g\n Chosen: %s\n', value, num2str(chosen));
  end
  
  if (minvalue<opt.stopvalue)
    if opt.display
      fprintf('Objective function value smaller than specified goal\n');
    end
    break
  end
end
x=rec(minvaluei).chosen;

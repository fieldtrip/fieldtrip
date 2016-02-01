function [x,rec] = bsearch(fun,x0,opt,varargin)
% BSEARCH Finds the minimum of a combinatorial function using backward search
%
%   X = BSEARCH(FUN, X0) attempts to find a combination of elements
%   of X0 which locally minimize the function FUN using backward
%   search (backward elimination). FUN accepts input X and returns
%   scalar function value F evaluated at X. X0 must be a vector of
%   initial indexes. Returned X contains indexes locally minimizing
%   the function FUN.
%
%   X = BSEARCH(FUN, X0, OPTIONS) allows use of optional
%   search parameters. See BSEARCH_OPT for details.
%
%   X = BSEARCH(FUN, X0, OPTIONS, P1, ..., Pn) P1,...,Pn are
%   additional parameters passed to the function FUN.
%
%   [X,REC] = BSEARCH(FUN, X0, ...) returns record of search as a
%   array of structs.
%
%   In backward search the elements are removed one at the time. 
%   The element ehich removal minimizes the function value is
%   removed at each round. In a case of tie choose randomly.
%

% Copyright (c) Aki Vehtari (2007)

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% 2007-02-23  Aki Vehtari  <Aki.Vehtari@hut.fi>
%             Use same parameter order as in fminunc etc.
%             Minimize instead of maximize
%             Rewrite.

% Columns
nx = size(x0,2);  % number of elements
nremoved = 0;     % number of removed elements
% Options
opt = bsearch_opt(opt);

% Base value
value = feval(fun,x0,varargin{:});
if (opt.display >= 1)
  fprintf(' Base value: %.4g\n Elements: %s\n', value, num2str(x0));
end
minvalue=value;
rec(nremoved+1).chosen = x0;
rec(nremoved+1).candidates = [];
rec(nremoved+1).values = value;

% The loop     
while (nx > opt.nsel)
  values=zeros(nx,1);
  for i1=1:nx
    values(i1) = feval(fun,setdiff(x0,x0(i1)),varargin{:});
    if (opt.display >= 2)
      fprintf('%d %.5g\n', x0(i1),  values(i1));
    end
  end
  % Let's get the index of the best column
  % In a case of tie, choose randomly
  value=min(values);
  mini=randpick(find(values==value));
  if (opt.stop && value>=minvalue)
    if opt.display
      fprintf('Stopping because value not decreasing\n');
    end
    break
  end
  nremoved = nremoved + 1;
  rec(nremoved+1).candidates = x0;
  rec(nremoved+1).values = values;
  x0 = setdiff(x0,x0(mini));       % Remove the element from x0
  nx=nx-1;
  rec(nremoved+1).chosen = x0;
  if value<minvalue
    minvalue=value;
    minvaluei=nremoved+1;
  end

  if (opt.display >= 1)
    fprintf(' Value: %.4g\n Chosen: %s\n', value, num2str(x0));
  end
  
  if (minvalue<opt.stopvalue)
    if opt.display
      fprintf('Objective function value smaller than specified goal\n');
    end
    break
  end
end
x=rec(minvaluei).chosen;

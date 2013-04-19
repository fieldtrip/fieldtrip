function s=binsgeq(pc,p)
%BINSGEQ  Binary search of sorted vector
%
%  Description:
%    S=BINSGEQ(PC,P) Binary search of a sorted (in ascending order) 
%    vector PC for a first element which is greater than or equal to 
%    a key P. If no greater element is found, returns LENGTH(PC). 
%    The size of S is the size of P.

% Copyright (c) 2004 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

% This is fast only in Matlab 6.5 and later with JIT accelerator
n=numel(pc);
[pm,pn]=size(p);
s=zeros(pm,pn);
for i=1:(pm*pn)
  pi=p(i);
  if pi<=pc(1)
    s(i)=1;
  else
    low=1;
    high=n;
    while (high-low)>1
      mid=floor((high+low)/2);
      if pi>pc(mid)
        low=mid;
      else
        high=mid;
      end
    end
    s(i)=high;
  end
end

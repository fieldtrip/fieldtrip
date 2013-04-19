function p1 = gibbs(p1, x)
%GIBBS      Gibbs sampling
%
%           Description
%           P1 = GIBBS(P1, X) Samples one sample for parameters
%           contained in structure P1 at level n and and gets a 
%           parameter structure X from level n-1 in hierarchy.
%           Returns a structure P1 containing new samples for
%           parameters in structure P1. 
%
%           If P1 contains higher level parameters, e.g P1 =
%           net.p.w{1} and though there is a field P1.p, the
%           function calls itself to sample for those higher level
%           parameters also. So function takes care of the whole
%           hierarchical structure of ARD if called for the lowest
%           level parameters. 

% Copyright (c) 1999-2003 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.


%Matlab 6.5
%if any(strcmp(fieldnames(p1),'p')) & ~isempty(p1.p) %isfield(p1,'p')
%  for a = fieldnames(p1.a)';
%    a=a{:};
%    if any(strcmp(fieldnames(p1.p),a)) %isfield(p1.p,a);
%      p2=p1.p.(a);
%      % call the spesified distribution sampler
%      p1.a.(a)=feval(['cond_' p1.f '_' p2.f],a,p1.a,p2.a,x);
%      % sample the upper level parameters
%      p1.p.(a)=gibbs(p2, p1.a.(a));
%    end
%  end
%end

% Pre matlab 6.5
S.type='.';
if any(strcmp(fieldnames(p1),'p')) & ~isempty(p1.p) %isfield(p1,'p')
  for a = fieldnames(p1.a)';
    a=a{:};
    S.subs=a;
    if any(strcmp(fieldnames(p1.p),a)) %isfield(p1.p,a);
      p2=subsref(p1.p,S);
      p1.a=subsasgn(p1.a, S, feval(['cond_' p1.f '_' p2.f],a,p1.a,p2.a,x)); 
      p1.p=subsasgn(p1.p, S, gibbs(p2, subsref(p1.a,S)));
    end
  end
end

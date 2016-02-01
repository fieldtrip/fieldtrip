function ea = ext_auc(P,tt,t)
%EXT_AUC Compute Extended AUC proposed by Chambless et al (2011)
%
%  Description
%    EA=EXT_AUC(P,TT,T) Compute Extended AUC(T) proposed by
%    Chambless et al (2011) given probability matrix P of size n x
%    size(tt,2), time vector TT and time T. Time T has to be one of
%    values in time vector TT.
%
%  Reference
%    L. E. Chambless, C. P. Cummiskey, and G. Cui (2011). Several
%    methods to assess improvement in risk prediction models:
%    Extension to survival analysis. Statistics in Medicine
%    30(1):22-38.
  
% Copyright (C) 2012 Ernesto Ulloa, Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.
  
  ip=inputParser;
  ip.addRequired('P',@(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('tt', @(x) isreal(x) && all(isfinite(x(:))))
  ip.addRequired('t', @(x) isreal(x) && all(isfinite(x(:))))
  ip.parse(P,tt,t);
  
  [n,nin]=size(P);
  S=1-P;
  D=-diff([; ones(n,1) S],1,2);
  sd=tt(2)-tt(1);
  
  if isempty(find(tt==t,1))
    error('Could not find your time value in the model prediction matrix')
  else
    I=find(tt==t);
  end
  
  if size(tt,2) > size(P,2)
    error('prediction matrix P must have number of columns equal to length(tt)')
  end
  
  comp=bsxfun(@gt,P(:,I),P(:,I)');
  %inner=D*S';
  den=(2/(1-mean(S(:,I)).^2))*(1/n^2);
  %ea=den*sum(sum(inner.*comp));
  extauc=0;
  for i=1:n
    for j=1:n
      if comp(i,j) == 1
        for d=1:size(tt,2)
          extauc=extauc + sum(D(i,d).*S(j,d)) + (sd/2)*(D(i,d));
        end
      end
    end
  end
  
  ea=extauc*den;

end


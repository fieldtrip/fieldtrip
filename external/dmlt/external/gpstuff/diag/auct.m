function at = auct(crit,y,z,tt)
%AUCS Compute area under curve for survival model at given time
%
%  Description
%    A = AUCT(CRIT,Z) Compute are under curve for survival model at
%    given time using criteria vector CRIT (where larger value
%    means larger risk of incidence), observed time vector Y,
%    censoring indicator matrix Z at times TT (0=event, 1=censored)
%    and a time vector TT.
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
ip.addRequired('crit',@(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) isreal(x) && all(isfinite(x(:))))
ip.addRequired('z', @(x) isreal(x) && all(isfinite(x(:))))
ip.addRequired('tt', @(x) isreal(x) && all(isfinite(x(:))))
ip.parse(crit,y,z,tt)

for i=1:size(tt,2)
  comp=bsxfun(@times,bsxfun(@and,y(:,i)<=tt(i),1-z(:,i)),bsxfun(@or,bsxfun(@and,y(:,i)<=tt(i),z(:,i)),y(:,i)>=tt(i))');
  conc=bsxfun(@times,bsxfun(@gt,crit(:,i),crit(:,i)'),comp);
  at(i,1)=sum(conc(:))./sum(comp(:));
end

end


function ct = hct(crit,y,z,tt)
%HCT Compute Harrel's C for survival model at several time points
%
%  Description
%    CT = HCT(CRIT,Y,Z,TT) Compute Harrel's C statistic estimate at
%    times TT using criteria vector CRIT (where larger value means
%    larger risk of incidence), observed time vector Y and event
%    indicator vector Z (0=event, 1=censored)
%
%  Reference
%    L. E. Chambless, C. P. Cummiskey, and G. Cui (2011). Several
%    methods to assess improvement in risk prediction models:
%    Extension to survival analysis. Statistics in Medicine
%    30(1):22-38.
%

% Copyright (C) 2012 Ernesto Ulloa, Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

ip=inputParser;
ip.addRequired('crit',@(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y',@(x) ~isempty(x) && isreal(x))
ip.addRequired('z', @(x) ~isempty(x) && isreal(x))
ip.addRequired('tt', @(x) ~isempty(x) && isreal(x))

ip.parse(crit,y,z,tt)

if size(y,2) ~= size(z,2)
  error('y and z dimensions must match')   
end

for i=1:size(tt,2)
  comp=bsxfun(@and,bsxfun(@and,bsxfun(@lt,y(:,i),y(:,i)'),y(:,i)<=tt(i)),z(:,end)==0);
  conc=bsxfun(@gt,crit(:,i),crit(:,i)').*comp;
  ct(i,1)=sum(conc(:))./sum(comp(:));
end

end


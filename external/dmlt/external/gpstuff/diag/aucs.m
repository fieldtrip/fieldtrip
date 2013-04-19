function [a,fps,tps] = aucs(crit,z)
%AUCS Compute area under curve for survival model
%
%  Description
%    A = AUCS(CRIT,Z) Compute are under curve for survival model
%    using criteria vector CRIT (where larger value means larger
%    risk of incidence) and censoring indicator column vector Z
%    (0=event, 1=censored).
%
%    [A,FPS,TPS] = AUCS(CRIT,Z) Return also a vector of false
%    positive rates FPS and a vector of true positive rates TPS.
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
ip.addRequired('z', @(x) isreal(x) && all(isfinite(x(:))))
ip.parse(crit,z)

ye=z;
tps=0;
fps=0;
qs=sort(crit,'descend');

for i2=1:numel(crit)
  tps(i2+1)=mean(crit>=qs(i2) & 1-ye)./mean(1-ye);
  fps(i2+1)=mean(crit>=qs(i2) & ye)./(mean(ye));
end

if (mean(1-ye)==0)||(mean(ye)==0) 
  warning('z vector has no different values, function will return 0');
  a=0;
else
  a=sum([diff(fps).*mean([tps(1:end-1);tps(2:end)],1)]);  
end

end

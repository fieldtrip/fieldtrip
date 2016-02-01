function [au,c] = assess(crit1,crit2,y,z,tt)
%   [AU,C] = ASSESS(CRIT1,CRIT2,Y,Z,TT)
%
%   Given two criteria vectors for two different models (CRIT1 and CRIT2), 
%   observed time matrix Y, event indicator matrix z (= 0 if event is experienced before tt and  =1 if not)
%   and time vector tt returns Harrel's C and AUC for each value in TT time vector
%   
%   Note: Z dimensions must be size(y,1) X size(tt,2)  

% Copyright (C) 2012 Ernesto Ulloa, Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

ip=inputParser;
ip.addRequired('crit1',@(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('crit2',@(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) isreal(x) && all(isfinite(x(:))))
ip.addRequired('z', @(x) isreal(x) && all(isfinite(x(:))))
ip.addRequired('tt', @(x) isreal(x) && all(isfinite(x(:))))
ip.parse(crit1,crit2,y,z,tt)
 
    for i=1:size(tt,2)
        c1(i)=hct(crit1(:,i),y(:,i),z(:,size(tt,2)),tt(i));
        c2(i)=hct(crit2(:,i),y(:,i),z(:,size(tt,2)),tt(i));
        au1(i)=auct(crit1(:,i),y(:,i),z(:,i),tt(i));
        au2(i)=auct(crit2(:,i),y(:,i),z(:,i),tt(i));
    end
    
   c(:,1)=c1;
   c(:,2)=c2;
   au(:,1)=au1;
   au(:,2)=au2;
end


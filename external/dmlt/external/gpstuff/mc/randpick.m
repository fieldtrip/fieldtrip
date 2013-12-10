function r = randpick(x)
% RANDPICK - Pick element from x randomly
%            If x is matrix, pick row from x randomly.

%   Author: Aki Vehtari <Aki.Vehtari@hut.fi>
%   Last modified: 2004-09-07 11:25:23 EEST

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

if any(size(x)==1)
  r=x(floor(rand.*length(x)+1));
else
  r=x(floor(rand.*size(x,1)+1),:);
end

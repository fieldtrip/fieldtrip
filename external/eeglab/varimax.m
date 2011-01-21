% varimax() - Perform orthogonal Varimax rotation on rows of a data 
%             matrix.
%
% Usage: >> V = varimax(data); 
%        >> [V,rotdata] = varimax(data,tol);
%        >> [V,rotdata] = varimax(data,tol,'noreorder') 
%
% Inputs:
%   data        - data matrix
%   tol         - set the termination tolerance to tol {default: 1e-4}
%   'noreorder' - Perform the rotation without component reorientation 
%                 or reordering by size. This suppression is desirable 
%                 when doing a q-mode analysis. {default|0|[] -> reorder}
% Outputs:
%   V           - orthogonal rotation matrix, hence
%   rotdata     - rotated matrix, rotdata = V*data;
%
% Author: Sigurd Enghoff - CNL / Salk Institute, La Jolla 6/18/98
%
% See also: runica(), pcasvd(), promax()

% Copyright (C) Sigurd Enghoff - CNL / Salk Institute, La Jolla 6/18/98
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% Reference: % Henry F. Kaiser (1958) The Varimx criterion for 
% analytic rotation in factor analysis. Pychometrika 23:187-200.
%
% modified to return V alone by Scott Makeig, 6/23/98
% 01-25-02 reformated help & license, added link -ad 

function [V,data] = varimax(data,tol,reorder)

if nargin < 1
   help varimax
   return
end

DEFAULT_TOL = 1e-4;  % default tolerance, for use in stopping the iteration
DEFAULT_REORDER = 1; % default to reordering the output rows by size
                     % and adjusting their sign to be rms positive.
MAX_ITERATIONS = 50; % Default
qrtr = .25;          % fixed value

if nargin < 3
    reorder = DEFAULT_REORDER;
elseif isempty(reorder) | reorder == 0
    reorder = 1; % set default
else
    reorder = strcmp('reorder',reorder);
end

if nargin < 2
    tol = 0;
end
if tol == 0
   tol = DEFAULT_TOL;
end
if ischar(tol)
   fprintf('varimax(): tol must be a number > 0\n');
   help varimax
   return
end

eps1 = tol; % varimax toler
eps2 = tol;

V = eye(size(data,1)); % do unto 'V' what is done to data
crit = [sum(sum(data'.^4)-sum(data'.^2).^2/size(data,2)) 0];
inoim = 0;
iflip = 1;
ict = 0;

fprintf(...
  'Finding the orthogonal Varimax rotation using delta tolerance %d...\n',...
                                                                eps1);
while inoim < 2 & ict < MAX_ITERATIONS & iflip,
    iflip = 0;
    for j = 1:size(data,1)-1,
        for k = j+1:size(data,1),
            u = data(j,:).^2-data(k,:).^2;
            v = 2*data(j,:).*data(k,:);
            a = sum(u);
            b = sum(v);
            c = sum(u.^2-v.^2);
            d = sum(u.*v);

            fden = size(data,2)*c + b^2 - a^2;
            fnum = 2 * (size(data,2)*d - a*b);

            if abs(fnum) > eps1*abs(fden)
                iflip = 1;
                angl = qrtr*atan2(fnum,fden);
                tmp    =  cos(angl)*V(j,:)+sin(angl)*V(k,:);
                V(k,:) = -sin(angl)*V(j,:)+cos(angl)*V(k,:);
                V(j,:) = tmp;

                tmp       =  cos(angl)*data(j,:)+sin(angl)*data(k,:);
                data(k,:) = -sin(angl)*data(j,:)+cos(angl)*data(k,:);
                data(j,:) = tmp;
            end
        end
    end

    crit = [sum(sum(data'.^4)-sum(data'.^2).^2/size(data,2)) crit(1)];
    inoim = inoim + 1;
    ict = ict + 1;

    fprintf('#%d - delta = %g\n',ict,(crit(1)-crit(2))/crit(1));

    if (crit(1) - crit(2)) / crit(1) > eps2
        inoim = 0;
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if reorder
    fprintf('Reordering rows...');
    [fnorm index] = sort(sum(data'.^2));
    V = V .* ((2 * (sum(data') > 0) - 1)' * ones(1, size(V,2)));
    data = data .* ((2 * (sum(data') > 0) - 1)' * ones(1, size(data,2)));
    V = V(fliplr(index),:);
    data = data(fliplr(index),:);
    fprintf('\n');
else
    fprintf('Not reordering rows.\n');
end

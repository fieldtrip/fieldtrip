% nanstd() - std, not considering NaN values
%
% Usage: same as std()

% Author: Arnaud Delorme, CNL / Salk Institute, Sept 2003

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% $Log: nanstd.m,v $
% Revision 1.2  2006/08/19 12:27:02  marsie
% fixed bug: crashed with empty flag
%
% Revision 1.1  2006/03/20 14:36:46  jansch
% adjusted the nan_XXX functions in fieldtrip's private directory such that the
% corresponding functions behave consistently with the identical matlab-functions,
% not using the stats toolbox. original private functions are renamed from nan_XXX
% into nanXXX (also consistent with matlab terminology)
%
% Revision 1.1  2005/04/08 06:53:00  roboos
% originates from EEGLAB, goes together with nan_mean
%
% Revision 1.1  2003/09/04 00:57:11  arno
% Initial revision
%

function out = nanstd(in, varargin)
   
if nargin < 1
  help nanstd;
  return;
end
if nargin == 1, flag = 0; end
if nargin <  3,
  if size(in,1) ~= 1
    dim = 1;
  elseif size(in,2) ~= 1 
    dim = 2; 
  else
    dim = 3;
  end
end
if nargin == 2, flag = varargin{1}; end
if nargin == 3,
  flag = varargin{1};
  dim  = varargin{2};
end
if isempty(flag), flag = 0; end

nans = find(isnan(in));
in(nans) = 0;
   
nonnans = ones(size(in));
nonnans(nans) = 0;
nonnans = sum(nonnans, dim);
nononnans = find(nonnans==0);
nonnans(nononnans) = 1;
   
out = sqrt((sum(in.^2, dim)-sum(in, dim).^2./nonnans)./(nonnans-abs(flag-1)));
out(nononnans) = NaN;

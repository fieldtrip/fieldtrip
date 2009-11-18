function s = mynanstd(data)
%MYNANSTD take standard deviation along first dimensions ignoring nans
%
% Author: Marcel van Gerven, based on:
%
% Author: Arnaud Delorme, CNL / Salk Institute, Sept 2003
%
%123456789012345678901234567890123456789012345678901234567890123456789012
%
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
%
%   $Log: mynanstd.m,v $
%
    nans = (isnan(data));
    data(nans) = 0;
   
    nonnans = ones(size(data));
    nonnans(nans) = 0;
    nonnans = sum(nonnans);
    nononnans = (nonnans==0);
    nonnans(nononnans) = 1;
   
    s = sqrt((sum(data.^2)-sum(data).^2./nonnans)./(nonnans-1));
    s(isnan(s)) = 0; % arises if we have only one data point
    
    s(nononnans) = nan;   

end

function [num,dims]=dimnum(dimord, dim)
% This function returns the number of the given dimention 'dim' in the dimord string.
%
% Syntax: [num,dims]=dimnum(dimord, dim)
%
% e.g. when dimord='rpt_chancmb_freq_time' and dim='time', dimnum returns num=4
%       and dims contains {'rpt','chancmb','freq','tim'}.
% e.g. when dimord='rpt_chancmb_freq_time' and dim='chancmb', dimnum returns num=2
%       and dims again contains {'rpt','chancmb','freq','tim'}.
% 
% For the known dimentiontypes dim can also be 'time' or 'frequency'.
% The known types are:
% tim: 'time'
% freq: 'frq', 'frequency'
% chancmb: 'sgncmb', 'channel', 'signal combination', 'channels'
% rpt: 'trial','trials'
%
% When dim is not found in dimord, an empty matrix is returned, but
% dims then still contains all dims in dimord.

% Copyright (C) 2005, Geerten Kramer
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if(~ischar(dimord)||~ischar(dim))error('Both dimord and dim must be strings');end;

dims=tokenize(dimord,'_'); % splits the dimord string in parts using '_' as a delimiter
dim=lower(dim); % makes the function case unsensitive.

switch dim %convert some dimension names.
case 'tim'
    dim='time';
case {'frq','frequency'}
    dim='freq';
case {'channel','signal combination','channels','sgncmb','chan'}
    %dim='chancmb';
    dim='chan';
case {'trial','trials'}
    dim='rpt';
end

num=find(strcmp(dims,dim)); % find the number of the specified dimension.




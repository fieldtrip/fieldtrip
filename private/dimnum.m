function [num,dims]=dimnum(dimord, dim)
% This function returns the number of the given dimention 'dim' in the dimord string.
%
% Syntax: [num,dims]=dimnum(dimord, dim)
%
% e.g. when dimord='rpt_chancmb_freq_time' and dim='time', dimnum returns num=4
% 		and dims contains {'rpt','chancmb','freq','tim'}.
% e.g. when dimord='rpt_chancmb_freq_time' and dim='chancmb', dimnum returns num=2
% 		and dims again contains {'rpt','chancmb','freq','tim'}.
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
% $Log: dimnum.m,v $
% Revision 1.4  2009/02/02 13:21:51  jansch
% changed 'chancmb' into 'chan', added 'chan'
%
% Revision 1.3  2006/02/01 12:26:04  roboos
% made all uses of dimord consistent with the common definition of data dimensions, see the fixdimord() function
%
% Revision 1.2  2005/10/04 12:19:26  geekra
% Added comment and discription.
%
% Revision 1.1  2005/10/04 12:06:00  geekra
% Initial implementation.
%

if(~isstr(dimord)||~isstr(dim))error('Both dimord and dim must be strings');end;

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
end;

num=find(strcmp(dims,dim)); % find the number of the specified dimension.



